//
//  main.cpp
//  FEXIPRO
//
//  Created by Firas Abuzaid on 3/17/17.
//  Copyright © 2017 FutureData. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <armadillo>

namespace chrono = std::chrono;
namespace opt = boost::program_options;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;

double *parse_weights_csv(const std::string filename, const int num_rows,
                         const int num_cols) {
  std::cout << "Loading " << filename << "...." << std::endl;
  std::ifstream weight_file(filename.c_str(), std::ios_base::in);
  double *weights = new double[num_rows * num_cols];

  std::string buffer;
  if (weight_file) {
    for (int i = 0; i < num_rows; i++) {
      double *d = &weights[i * num_cols];
      for (int j = 0; j < num_cols; j++) {
        double f;
        weight_file >> f;
        if (j != num_cols - 1) {
          std::getline(weight_file, buffer, ',');
        }
        d[j] = f;
      }
      std::getline(weight_file, buffer);
    }
  }
  weight_file.close();
  return weights;
}

arma::mat parse_weights_csv(const std::string filename) {
  std::cout << "Loading " << filename << "...." << std::endl;
  arma::mat weights;
  const bool success = weights.load(filename, arma::csv_ascii);
  if (success) {
    std::cout << filename << " successfully loaded" << std::endl;
  }
  return weights.t();
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
      "Nubmer of latent factors")/*(
      "num-clusters,c", opt::value<size_t>()->required(), "Number of clusters")(
      "sample-ratio,s", opt::value<float>()->default_value(0.1),
      "Ratio of users to sample during clustering, between 0. and 1.")(
      "num-iters,i", opt::value<size_t>()->default_value(3),
      "Number of iterations to run clustering, default: 10")
      ("num-bins,b", opt::value<size_t>()->default_value(1), "Number of bins,
      default: 1")
      ("num-threads,t", opt::value<size_t>()->default_value(1),
       "Number of threads, default: 1")*/;

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
  // const size_t num_clusters = args["num-clusters"].as<size_t>();
  // const float sample_ratio = args["sample-ratio"].as<float>();
  // const size_t num_iters = args["num-iters"].as<size_t>();
  // const size_t num_bins = args["num-bins"].as<size_t>();
  // const size_t num_threads = args["num-threads"].as<size_t>();

  // 1) Load user and item weights
  arma::mat user_weights = parse_weights_csv(user_weights_file);
  std::cout << "User matrix: " << arma::size(user_weights) << std::endl;

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

  arma::mat item_weights = parse_weights_csv(item_weights_file);
  std::cout << "Item matrix: " << arma::size(item_weights) << std::endl;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  auto t2 = Time::now();
  const bool success = svd_econ(U, s, V, item_weights, "both", "dc");
  auto t3 = Time::now();

  fsec thin_svd_time_s = t3 - t2;

  ms thin_svd_time_ms = std::chrono::duration_cast<ms>(thin_svd_time_s);

  std::cout << "U: " << arma::size(U) << std::endl;
  std::cout << "s: " << arma::size(s) << std::endl;
  std::cout << "V: " << arma::size(V) << std::endl;
  std::cout << success << std::endl;

  std::cout << "Thin SVD:" << std::endl;
  std::cout << thin_svd_time_s.count() << "s\n";
  std::cout << thin_svd_time_ms.count() << "ms\n";

  return 0;
}
