//
//  main.cpp
//  FEXIPRO
//
//  Created by Firas Abuzaid on 3/17/17.
//  Copyright © 2017 FutureData. All rights reserved.
//

#include <armadillo>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <omp.h>

namespace chrono = std::chrono;
namespace opt = boost::program_options;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;

const double RHO = 0.7;
const int E = 100;

// For debugging purposes in gdb, since we can't use overloaded operators
double get_elem(arma::mat const &m, int i, int j) { return m(i, j); }
double get_elem(arma::vec const &v, int i) { return v(i); }
int get_elem(arma::uvec const &v, int i) { return v(i); }

arma::mat parse_weights_csv(const std::string filename) {
  std::cout << "Loading " << filename << "...." << std::endl;
  arma::mat weights;
  const bool success = weights.load(filename, arma::csv_ascii);
  if (success) {
    std::cout << filename << " successfully loaded" << std::endl;
  }
  return weights.t();
}

// Preprocess items matrix (P) according to Algorithm 3 in the FEXIPRO paper
void preprocess(const int d, arma::mat &P, int &w, arma::mat &U, arma::vec &s,
                arma::vec &p_norms, arma::vec &p_bar_h_norms,
                arma::vec &p_double_hat_h_norms) {
  p_norms = arma::zeros(P.n_cols);
  p_bar_h_norms = arma::zeros(P.n_cols);
  p_double_hat_h_norms = arma::zeros(P.n_cols);

  // 1) sort item weights by decreasing length
  for (int i = 0; i < P.n_cols; ++i) {
    p_norms(i) = arma::norm(P.col(i), 2);
  }
  arma::uvec indices = arma::sort_index(p_norms, "descend");
  P = P.cols(indices);

  // 2) compute thin SVD to get U, Sigma_d, V_1
  arma::mat V;
  const bool success = svd_econ(U, s, V, P, "both", "dc");

  // 3) Calculate w according to rho
  const double total = arma::sum(s);
  double sum = 0.0;
  w = 0;
  for (; w < s.n_elem; ++w) {
    // std::cout << s(w) << std::endl;
    sum += s(w);
    if (sum / total > RHO) {
      break;
    }
  }
  // std::cout << w << std::endl;
  // 4)
  V = V.t();
  const double max_P_l = P.head_rows(w).max();
  const double max_P_h = P.tail_rows(d - w).max();

  const double b = p_norms.max();
  const double c_s_init = std::max(1.0, std::abs(P.min()));
  arma::vec c = arma::zeros(d);
  for (int i = 0; i < d; ++i) {
    c(i) = c_s_init + s(i) / s(d - 1);
  }

  for (int i = 0; i < V.n_cols; ++i) {
    const arma::vec p = P.col(i);
    const arma::vec p_bar = V.col(i);
    arma::vec p_hat_l = p_bar.head(w);
    arma::vec p_hat_h = p_bar.tail(d - w);
    p_bar_h_norms(i) = arma::norm(p_hat_h, 2);
    p_hat_l /= max_P_l;
    p_hat_h /= max_P_h;
    p_hat_l *= E;
    p_hat_h *= E;

    arma::vec p_hat_floor = arma::floor(p_bar);
    const double p_norm = p_norms(i);

    arma::vec p_double_hat = arma::zeros(d + 2);
    arma::vec p_tilde = arma::zeros(d + 1);
    p_tilde(0) = std::sqrt(b * b - p_norm * p_norm);
    p_tilde.tail(d) = p + c;
    p_double_hat(0) = arma::norm(p_tilde, 2);
    p_double_hat.tail(d + 1) = p_tilde;
    p_double_hat_h_norms(i) = arma::norm(p_double_hat.tail(d - w), 2);
  }
}

int main(int argc, const char *argv[]) {
  opt::options_description description("SimDex");
  description.add_options()("help,h", "Show help")(
      "weights-dir,w", opt::value<std::string>()->required(),
      "weights directory; must contain user_weights.csv and item_weights.csv")(
      "top-k,k", opt::value<size_t>()->required(),
      "Top K items to return per user")(
      "num-users,m", opt::value<size_t>()->required(), "Number of users")(
      "num-items,n", opt::value<size_t>()->required(), "Number of items")(
      "num-latent-factors,f", opt::value<size_t>()->required(),
      "Number of latent factors");

  opt::variables_map args;
  opt::store(opt::command_line_parser(argc, argv).options(description).run(),
             args);

  if (args.count("help")) {
    std::cout << description << std::endl;
    exit(1);
  }

  opt::notify(args);

  const std::string weights_dir = args["weights-dir"].as<std::string>();
  const std::string user_weights_file = weights_dir + "/user_weights.csv";
  const std::string item_weights_file = weights_dir + "/item_weights.csv";

  const size_t K = args["top-k"].as<size_t>();
  const size_t num_users = args["num-users"].as<size_t>();
  const size_t num_items = args["num-items"].as<size_t>();
  const size_t num_latent_factors = args["num-latent-factors"].as<size_t>();
  omp_set_num_threads(1);  // always set to one thread

  // Load item weights
  arma::mat item_weights = parse_weights_csv(item_weights_file);
  std::cout << "Item matrix: " << arma::size(item_weights) << std::endl;

  // Preprocessing return values
  int w;
  arma::mat U;
  arma::vec s;
  arma::vec p_norms;
  arma::vec p_bar_h_norms;
  arma::vec p_double_hat_h_norms;
  // Preprocess
  auto t2 = Time::now();
  preprocess(num_latent_factors, item_weights, w, U, s, p_norms, p_bar_h_norms,
             p_double_hat_h_norms);
  auto t3 = Time::now();

  fsec preprocess_time_s = t3 - t2;
  ms preprocess_time_ms = std::chrono::duration_cast<ms>(preprocess_time_s);
  std::cout << "Preprocess time" << std::endl;
  std::cout << preprocess_time_s.count() << "s\n";
  std::cout << preprocess_time_ms.count() << "ms\n";

  // Load user weights
  arma::mat user_weights = parse_weights_csv(user_weights_file);
  std::cout << "User matrix: " << arma::size(user_weights) << std::endl;

  return 0;
}

// Clustering code for SimDex
//
// auto t0 = Time::now();
// arma::mat means;
// arma::kmeans(means, user_weights.submat(0, 0, num_latent_factors - 1,
//                                         int(num_users * sample_ratio)),
//              num_clusters, arma::static_subset, num_iters, false);
// auto t1 = Time::now();
// fsec clustering_time_s = t1 - t0;
// ms clustering_time_ms = std::chrono::duration_cast<ms>(clustering_time_s);
// std::cout << "Clustering:" << std::endl;
// std::cout << clustering_time_s.count() << "s\n";
// std::cout << clustering_time_ms.count() << "ms\n";
