//
//  main.cpp
//  FEXIPRO
//
//  Created by Firas Abuzaid on 3/17/17.
//  Copyright © 2017 FutureData. All rights reserved.
//

#include <armadillo>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#ifdef DEBUG
#include <cassert>
#endif
#include <chrono>
#include <float.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <queue>

namespace chrono = std::chrono;
namespace opt = boost::program_options;

using namespace arma;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;

const double RHO = 0.7;
const int E = 100;

static uint32_t full_dot_products = 0;

// For debugging purposes in gdb, since we can't use overloaded operators
#ifdef DEBUG
void print(const vec &v) { v.print(); }
void print(const ivec &v) { v.print(); }
double get_elem(const mat &m, int i, int j) { return m(i, j); }
double get_elem(const vec &v, int i) { return v(i); }
int get_elem(const uvec &v, int i) { return v(i); }
double get_norm(const vec &v) { return norm(v, 2); }
double get_dot(const vec &u, const vec &v) { return dot(u, v); }
#endif

mat parse_weights_csv(const std::string filename, const uint32_t scale) {
  std::cout << "Loading " << filename << "...." << std::endl;
  mat weights;
  const bool success = weights.load(filename, csv_ascii);
  if (success) {
    std::cout << filename << " successfully loaded" << std::endl;
  }
  return weights.t() / scale;
}

// implements Theorem 2 from FEXIPRO paper: IU(q, p) =
// \sum_{s=1}^d (floor(q_s)*floor(p_s) + abs(floor(q_s)) + abs(floor(p_s)) + 1)
// WARNING: args q and p should have already been floored
inline double integer_upper_bound(const ivec &q, const ivec &p,
                                  const uint32_t d) {
  return dot(q, p) + sum(abs(q)) + sum(abs(p)) + d;
}

#ifdef DEBUG
void check_double_hat_equivalency(const vec &q_double_hat,
                                  const vec &p_double_hat, const vec &q,
                                  const vec &p, const vec &c,
                                  const double q_norm,
                                  const double p_prime_squared_norm) {
  const double lhs = dot(q_double_hat, p_double_hat);
  const double rhs = 2 * dot(q, p) / q_norm +
                     2 * (dot(c, q) / q_norm + dot(c, p) + dot(c, c)) -
                     p_prime_squared_norm;
  const double epsilon = 1e-10;
  const double diff = lhs - rhs;
  assert((diff < epsilon) && (-diff < epsilon));
}
void check_double_hat_equivalency(const vec &q_double_hat,
                                  const vec &p_double_hat, const vec &q,
                                  const vec &p, const double q_norm,
                                  const double thresh_prime_online,
                                  const double thresh_prime_offline) {
  const double lhs = dot(q_double_hat, p_double_hat);
  const double rhs =
      2 * dot(q, p) / q_norm + thresh_prime_online + thresh_prime_offline;
  const double epsilon = 1e-10;
  const double diff = lhs - rhs;
  assert((diff < epsilon) && (-diff < epsilon));
}
#endif

// Preprocess items matrix (P) according to Algorithm 3 in the FEXIPRO paper
void preprocess(const uint32_t d, uint32_t &w, mat &P, mat &U, mat &P_bar,
                imat &P_hat_floors, uvec &item_ids, vec &sigma, vec &c,
                vec &p_norms, vec &p_bar_h_norms, vec &p_double_hat_h_norms,
                vec &thresh_prime_offline, vec &thresh_prime_offline_l,
                double &max_P_bar_l, double &max_P_bar_h
#ifdef DEBUG
                ,
                mat &P_double_hats, vec &p_prime_squared_norms
#endif
                ) {

  P_hat_floors = zeros<imat>(size(P));
  p_norms = zeros(P.n_cols);
  vec p_bar_norms = zeros(P.n_cols);
  p_bar_h_norms = zeros(P.n_cols);
  p_double_hat_h_norms = zeros(P.n_cols);
  thresh_prime_offline = zeros(P.n_cols);
  thresh_prime_offline_l = zeros(P.n_cols);

#ifdef DEBUG
  P_double_hats = zeros(d + 2, P.n_cols);
  p_prime_squared_norms = zeros(P.n_cols);
#endif

  // 1) sort item weights by decreasing length
  for (uint32_t i = 0; i < P.n_cols; ++i) {
    p_norms(i) = norm(P.unsafe_col(i), 2);
  }
  // sort item norms from highest to lowest, and keep
  // track of the original item ids
  item_ids = sort_index(p_norms, "descend");
  P = P.cols(item_ids);
  p_norms = p_norms(item_ids);

  // 2) compute thin SVD to get U, Sigma_d, V_1
  mat V;
  svd_econ(U, sigma, V, P, "both", "dc");

  // 3) Calculate w according to rho
  const double total = sum(sigma);
  double sum = 0.0;
  w = 0;
  for (; w < sigma.n_elem; ++w) {
    sum += sigma(w);
    if (sum / total > RHO) {
      ++w; // include current element
      break;
    }
  }

  // 4) calculate statistics for each item vector
  P_bar = V.t();

  // Compute b, which is needed for p_double_hat
  // b = max(||p||) for p in P_bar
  double b = -DBL_MAX;
  for (uint32_t i = 0; i < P_bar.n_cols; ++i) {
    const double p_bar_norm = norm(P_bar.unsafe_col(i), 2);
    if (p_bar_norm > b) {
      b = p_bar_norm;
    }
    // save in p_bar_norms, so we can use it later
    // to compute p_double_hat
    p_bar_norms(i) = p_bar_norm;
  }
  // Compute c, which is needed for p_double_hat
  // c = (max(1, |p_min|) + sigma_1/sigma_d, ..., max(1, |p_min|) +
  // sigma_d/sigma_d)
  // |p_min| = absolute value of minimum value in P_bar
  const double c_init = std::max(1.0, std::abs(P_bar.min()));
  c = zeros(d);

  const double threshold = 1e-10;
  double denom = 0.;
  for (uint32_t i = d-1; i > 0; --i) {
    if (sigma(i) < threshold) {
      continue;
    }
    denom = sigma(i);
    break;
  }

  for (uint32_t i = 0; i < d; ++i) {
    c(i) = c_init + sigma(i) / denom;
  }
  const vec c_l = c.head(w);

  max_P_bar_l = abs(P_bar.head_rows(w)).max();
  max_P_bar_h = abs(P_bar.tail_rows(d - w)).max();
  // should we use P or P_bar? P_bar, since P_bar = V^t, and we want to take
  // advantage of the SVD operation before we do integer pruning
  for (uint32_t i = 0; i < P_bar.n_cols; ++i) {
    const vec p_bar = P_bar.unsafe_col(i);

    //  1. Compute p_hat_l, p_hat_h, and p_hat_floor **for p_bar**
    // p_hat_l = (e*p_1/max(P_l), e*p_2/max(P_l), ..., e*p_w/max(P_l))
    vec p_hat_l = p_bar.head(w);
    p_hat_l /= max_P_bar_l;
    p_hat_l *= E;
    // p_hat_h = (e*p_w+1/max(P_h), e*p_w+2/max(P_h), ..., e*p_d/max(P_h))
    vec p_hat_h = p_bar.tail(d - w);
    // Calculate ||p_bar_h|| now, since we're going to need it later on
    p_bar_h_norms(i) = norm(p_hat_h, 2);
    p_hat_h /= max_P_bar_h;
    p_hat_h *= E;
    // p_hat_floor = floor of concat(p_hat_l, p_hat_h)
    const ivec p_hat_floor =
        conv_to<ivec>::from(floor(join_vert(p_hat_l, p_hat_h)));
    P_hat_floors.unsafe_col(i) = p_hat_floor;

    // 2. Compute p_double_hat and thresh_prime_offline
    // p_double_hat = (||p_prime||^2, p_prime_1, ..., p_prime_d+1)
    vec p_double_hat = zeros(d + 2);
    // p_prime = (sqrt(b^2 - ||p||^2), p_1 + c_1, ..., p_d + c_d)
    vec p_prime = zeros(d + 1);
    const double p_bar_norm = p_bar_norms(i);
    p_prime(0) = std::sqrt(b * b - p_bar_norm * p_bar_norm);
    p_prime.tail(d) = p_bar + c;
    double p_prime_squared_norm = norm(p_prime, 2);
    p_prime_squared_norm *= p_prime_squared_norm;
    p_double_hat(0) = p_prime_squared_norm;
    p_double_hat.tail(d + 1) = p_prime;
#ifdef DEBUG
    p_prime_squared_norms(i) = p_prime_squared_norm;
    P_double_hats.unsafe_col(i) = p_double_hat;
    // p_double_hat should always have nonnegative values
    for (int i = 0; i < d + 2; ++i) {
      assert(p_double_hat(i) >= 0.0);
    }
#endif
    // 3. Compute ||p||, ||p_bar_h||, and ||p_double_hat_h||
    // ||p|| was completed when we sorted P
    // ||p_bar_h|| was completed as we computed p_hat_h
    // ||p_double_hat_h|| is the only one missing
    p_double_hat_h_norms(i) = norm(p_double_hat.tail(d - w), 2);

    // thresh_prime_offline = 2*sum(c_1*p_1 + c_1^2, ..., c_d*p_d + c_d^2) -
    // ||p_prime||^2
    thresh_prime_offline(i) =
        2 * (dot(c, p_bar) + dot(c, c)) - p_prime_squared_norm;
    const vec p_bar_l = p_bar.head(w);
    thresh_prime_offline_l(i) =
        2 * (dot(c_l, p_bar_l) + dot(c_l, c_l)) - p_prime_squared_norm;
  }
}

// Algorithm 5 from the FEXIPRO paper: implement all types of pruning
double coordinate_scan(
    const vec &p_bar, const vec &q_bar, const ivec &p_hat_l_floor,
    const ivec &q_hat_l_floor, const ivec &p_hat_h_floor,
    const ivec &q_hat_h_floor, const vec &c_l, const uint32_t w,
    const uint32_t d, const double thresh, const double thresh_prime,
    const double q_bar_norm, const double p_bar_h_norm,
    const double q_bar_h_norm, const double p_double_hat_h_norm,
    const double q_double_hat_h_norm,
    const double thresh_prime_offline_l_single, const double max_P_bar_l,
    const double max_q_bar_l, const double max_P_bar_h, const double max_q_bar_h
#ifdef DEBUG
    ,
    const vec &p, const vec &q, const vec &c, const vec &p_double_hat,
    const vec &q_double_hat, const double p_prime_squared_norm,
    const double q_norm
#endif
    ) {
#ifdef DEBUG
  const double true_v = dot(p, q);
  const bool exit_early = true_v < thresh;
  double diff = 0.0;
  const double epsilon = 1e-10;
#endif

  const double b_l = (integer_upper_bound(q_hat_l_floor, p_hat_l_floor, w) *
                      max_q_bar_l * max_P_bar_l) /
                     (E * E);
  const double ub_1 = p_bar_h_norm * q_bar_h_norm;
  if (b_l + ub_1 < thresh) {
#ifdef DEBUG
    assert(exit_early);
#endif
    return -DBL_MAX;
  }
  const double b_h = (integer_upper_bound(q_hat_h_floor, p_hat_h_floor, d - w) *
                      max_q_bar_h * max_P_bar_h) /
                     (E * E);
  if (b_l + b_h < thresh) {
#ifdef DEBUG
    assert(exit_early);
#endif
    return -DBL_MAX;
  }
  const vec p_bar_l = p_bar.head(w);
  const vec q_bar_l = q_bar.head(w);
  double v = dot(p_bar_l, q_bar_l);
#ifdef DEBUG
  // dot(q_bar_l, p_bar_l) should always be less than
  // IU(q_bar_l, p_bar_l)
  assert(v < b_l);
#endif
  if (v + ub_1 < thresh) {
#ifdef DEBUG
    assert(exit_early);
#endif
    return -DBL_MAX;
  }
  // dot(q_double_hat_l, p_double_hat_l) = 2*dot(q_l, p_l)/||q|| +
  // 2*\sum_{s=1}^w (c_s*q_s/||q|| + c_s*p_s + c_s^2) - ||p_prime||^2
  const double dot_q_p_double_hat_l = 2 * v / q_bar_norm +
                                      2 * dot(c_l, q_bar_l) / q_bar_norm +
                                      thresh_prime_offline_l_single;
  const double ub_2 = p_double_hat_h_norm * q_double_hat_h_norm;
#ifdef DEBUG
  const vec p_double_hat_l = p_double_hat.head(w + 2);
  const vec q_double_hat_l = q_double_hat.head(w + 2);
  const double _dot_q_p_double_hat_l = dot(p_double_hat_l, q_double_hat_l);
  diff = _dot_q_p_double_hat_l - dot_q_p_double_hat_l;
  assert((diff < epsilon) && (-diff < epsilon));
#endif

  if (dot_q_p_double_hat_l + ub_2 < thresh_prime) {
#ifdef DEBUG
    assert(exit_early);
#endif
    return -DBL_MAX;
  }
  v += dot(q_bar.tail(d - w), p_bar.tail(d - w));
#ifdef DEBUG
  diff = v - true_v;
  assert((diff < epsilon) && (-diff < epsilon));
#endif
  full_dot_products++;
  return v;
}

// retrieve top K items for user q, based on Algorithm 4 in FEXIPRO paper
uvec retrieve_top_K(const uint32_t d, const uint32_t w, const uint32_t K,
                    const vec &q, const mat &P, const mat &P_bar, const mat &U,
                    const imat &p_hat_l_floors, const imat &p_hat_h_floors,
                    const uvec &item_ids, const vec &sigma, const vec &c,
                    const vec &c_l, const vec &p_norms,
                    const vec &p_bar_h_norms, const vec &p_double_hat_h_norms,
                    const vec &thresh_prime_offline,
                    const vec &thresh_prime_offline_l, const double max_P_bar_l,
                    const double max_P_bar_h
#ifdef DEBUG
                    ,
                    const mat &P_double_hats, const vec &p_prime_squared_norms
#endif
                    ) {

  uvec top_K_items = zeros<uvec>(K);

  std::priority_queue<std::pair<double, uint32_t>,
                      std::vector<std::pair<double, uint32_t> >,
                      std::greater<std::pair<double, uint32_t> > >
      queue;
  double thresh = -DBL_MAX;
  double thresh_prime = -DBL_MAX;

  // 1) Compute q_bar
  const vec q_bar = diagmat(sigma) * U.t() * q;
  // 2) Compute q_hat_l, q_hat_h, q_double_hat, and q_hat_floor

  // First, we compute q_hat_l and q_hat_h
  vec q_hat_l = q_bar.head(w);
  vec q_hat_h = q_bar.tail(d - w);
  const double max_q_bar_l = abs(q_hat_l).max();
  const double max_q_bar_h = abs(q_hat_h).max();
  q_hat_l /= max_q_bar_l;
  q_hat_h /= max_q_bar_h;
  q_hat_l *= E;
  q_hat_h *= E;

  // Then, we compute q_hat_l_floor and q_hat_h_floor
  const ivec q_hat_l_floor = conv_to<ivec>::from(floor(q_hat_l));
  const ivec q_hat_h_floor = conv_to<ivec>::from(floor(q_hat_h));

  // Finally, compute q_double_hat
  // q_double_hat = (-1, 2*q_prime_1, 2*q_prime_2, ..., 2*q_prime_d+1)
  vec q_double_hat = zeros(d + 2);
  q_double_hat(0) = -1;
  // q_prime = (0, q_1/||q|| + c_1, q_2/||q|| + c_2, ..., q_d/||q|| + c_d)
  vec q_prime = zeros(d + 1);
  q_prime.tail(d) = q_bar;
  q_prime = normalise(q_prime, 2);
  q_prime.tail(d) += c;
  q_double_hat.tail(d + 1) = 2 * q_prime;

#ifdef DEBUG
  // q_double_hat's first element should == -1; the rest should be nonnegative
  assert(q_double_hat(0) == -1);
  assert(q_double_hat(1) == 0);
  for (int i = 2; i < d + 2; ++i) {
    assert(q_double_hat(i) >= 0.0);
  }
#endif

  // 3) Compute ||q||, ||q_bar_h||, and ||q_double_hat_h||
  const double q_norm = norm(q, 2);
  const double q_bar_norm = norm(q_bar, 2);
  const double q_bar_h_norm = norm(q_bar.tail(d - w), 2);
  const double q_double_hat_h_norm = norm(q_double_hat.tail(d - w), 2);

  // Compute 2*dot(c, q)/||q||
  const double thresh_prime_online = 2 * dot(c, q_bar) / q_bar_norm;

  // 4) Search for top K
  uint32_t item_ind = 0;
  for (; item_ind < K; ++item_ind) {
    const vec p_bar = P_bar.unsafe_col(item_ind);
#ifdef DEBUG
    const vec p = P.unsafe_col(item_ind);
    const vec p_double_hat = P_double_hats.unsafe_col(item_ind);
    const double p_prime_squared_norm = p_prime_squared_norms(item_ind);
    const double thresh_prime_offline_single = thresh_prime_offline(item_ind);
    check_double_hat_equivalency(q_double_hat, p_double_hat, q_bar, p_bar, c,
                                 q_bar_norm, p_prime_squared_norm);
    check_double_hat_equivalency(
        q_double_hat.head(w + 2), p_double_hat.head(w + 2), q_bar.head(w),
        p_bar.head(w), c.head(w), q_bar_norm, p_prime_squared_norm);
    check_double_hat_equivalency(q_double_hat, p_double_hat, q_bar, p_bar,
                                 q_bar_norm, thresh_prime_online,
                                 thresh_prime_offline_single);
#endif
    const double score = dot(q_bar, p_bar);
    queue.push(std::make_pair(score, item_ind));
  }
  thresh = queue.top().first;
  thresh_prime = 2 * thresh / q_bar_norm + thresh_prime_online +
                 thresh_prime_offline(queue.top().second);

  for (item_ind = K; item_ind < P_bar.n_cols; ++item_ind) {
    const vec p_bar = P_bar.unsafe_col(item_ind);
    const vec p = P.unsafe_col(item_ind);
#ifdef DEBUG
    const vec p_double_hat = P_double_hats.unsafe_col(item_ind);
    const double p_prime_squared_norm = p_prime_squared_norms(item_ind);
    const double thresh_prime_offline_single = thresh_prime_offline(item_ind);
    check_double_hat_equivalency(q_double_hat, p_double_hat, q_bar, p_bar, c,
                                 q_bar_norm, p_prime_squared_norm);
    check_double_hat_equivalency(
        q_double_hat.head(w + 2), p_double_hat.head(w + 2), q_bar.head(w),
        p_bar.head(w), c.head(w), q_bar_norm, p_prime_squared_norm);
    check_double_hat_equivalency(q_double_hat, p_double_hat, q_bar, p_bar,
                                 q_bar_norm, thresh_prime_online,
                                 thresh_prime_offline_single);
#endif
    if (q_norm * p_norms(item_ind) <= thresh) {
      break;
    }
    const double v = coordinate_scan(
        p_bar, q_bar, p_hat_l_floors.unsafe_col(item_ind), q_hat_l_floor,
        p_hat_h_floors.unsafe_col(item_ind), q_hat_h_floor, c_l, w, d, thresh,
        thresh_prime, q_bar_norm, p_bar_h_norms(item_ind), q_bar_h_norm,
        p_double_hat_h_norms(item_ind), q_double_hat_h_norm,
        thresh_prime_offline_l(item_ind), max_P_bar_l, max_q_bar_l, max_P_bar_h,
        max_q_bar_h
#ifdef DEBUG
        ,
        p, q, c, P_double_hats.unsafe_col(item_ind), q_double_hat,
        p_prime_squared_norms(item_ind), q_norm
#endif
        );
    if (v > thresh) {
      queue.pop();
      queue.push(std::make_pair(v, item_ind));
      thresh = queue.top().first;
      thresh_prime = 2 * thresh / q_bar_norm + thresh_prime_online +
                     thresh_prime_offline(queue.top().second);
    }
  }
#ifdef DEBUG
  std::cout << "Full dot products: " << full_dot_products << std::endl;
#endif
  // 5) Get top K from priority queue
  for (uint32_t i = 0; i < K; i++) {
    const std::pair<double, int> pair = queue.top();
    top_K_items(K - i - 1) = item_ids(pair.second);
    queue.pop();
  }
  return top_K_items;
}

int main(int argc, const char *argv[]) {
  opt::options_description description("SimDex");
  description.add_options()("help,h", "Show help")(
      "weights-dir,w", opt::value<std::string>()->required(),
      "weights directory; must contain user_weights.csv and item_weights.csv")(
      "top-k,k", opt::value<uint32_t>()->required(),
      "Top K items to return per user")(
      "num-users,m", opt::value<uint32_t>()->required(), "Number of users")(
      "num-items,n", opt::value<uint32_t>()->required(), "Number of items")(
      "scale,s", opt::value<uint32_t>()->default_value(1), "scale weights by fraction")(
      "num-latent-factors,f", opt::value<uint32_t>()->required(),
      "Number of latent factors")("base-name",
                                  opt::value<std::string>()->required(),
                                  "Base name for output file");

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

  const uint32_t K = args["top-k"].as<uint32_t>();
  const uint32_t num_users = args["num-users"].as<uint32_t>();
  const uint32_t num_latent_factors = args["num-latent-factors"].as<uint32_t>();
  const uint32_t scale = args["scale"].as<uint32_t>();
  omp_set_num_threads(1);  // always set to one thread

  const std::string base_name = args["base-name"].as<std::string>();

  // Load item weights
  mat item_weights = parse_weights_csv(item_weights_file, scale);
  std::cout << "Item matrix: " << size(item_weights) << std::endl;

  // Preprocessing return values
  uint32_t w;
  mat U;                       // num_latent_factors x num_latent_factors
  mat P_bar;                   // num_latent_factors x num_items
  imat P_hat_floors;           // num_latent_factors x num_items
  uvec item_ids;               // num_items x 1
  vec sigma;                   // num_latent_factors x 1
  vec c;                       // num_latent_factors x 1
  vec p_norms;                 // num_items x 1
  vec p_bar_h_norms;           // num_items x 1
  vec p_double_hat_h_norms;    // num_items x 1
  vec thresh_prime_offline;    // num_items x 1
  vec thresh_prime_offline_l;  // num_items x 1
  double max_P_bar_l;
  double max_P_bar_h;
#ifdef DEBUG
  mat P_double_hats;
  vec p_prime_squared_norms;  // num_items x 1
#endif

  // Preprocess
  auto start = Time::now();
  preprocess(num_latent_factors, w, item_weights, U, P_bar, P_hat_floors,
             item_ids, sigma, c, p_norms, p_bar_h_norms, p_double_hat_h_norms,
             thresh_prime_offline, thresh_prime_offline_l, max_P_bar_l,
             max_P_bar_h
#ifdef DEBUG
             ,
             P_double_hats, p_prime_squared_norms
#endif
             );
  auto end = Time::now();

  fsec preprocess_time_s = end - start;
  ms preprocess_time_ms = std::chrono::duration_cast<ms>(preprocess_time_s);
  std::cout << "Preprocess time" << std::endl;
  std::cout << preprocess_time_s.count() << "s\n";
  std::cout << preprocess_time_ms.count() << "ms\n";

  // Load user weights
  mat user_weights = parse_weights_csv(user_weights_file, scale);
  std::cout << "User matrix: " << size(user_weights) << std::endl;
  start = Time::now();
  for (uint32_t i = 0; i < num_users; ++i) {
    const vec q = user_weights.unsafe_col(i);
    const uvec top_K_items = retrieve_top_K(
        num_latent_factors, w, K, q, item_weights, P_bar, U,
        P_hat_floors.head_rows(w),
        P_hat_floors.tail_rows(num_latent_factors - w), item_ids, sigma, c,
        c.head(w), p_norms, p_bar_h_norms, p_double_hat_h_norms,
        thresh_prime_offline, thresh_prime_offline_l, max_P_bar_l, max_P_bar_h
#ifdef DEBUG
        ,
        P_double_hats, p_prime_squared_norms
#endif
        );
#ifdef DEBUG
    std::cout << top_K_items(0) << std::endl;
#endif
  }
  end = Time::now();
  const double avg_full_dot_products = full_dot_products / (double)num_users;
  std::cout << "Avg. full dot products: " << avg_full_dot_products << std::endl;
  fsec comp_time_s = end - start;
  ms comp_time_ms = std::chrono::duration_cast<ms>(comp_time_s);
  std::cout << "Comp time" << std::endl;
  std::cout << comp_time_s.count() << "s\n";
  std::cout << comp_time_ms.count() << "ms\n";

  std::ofstream timing_stats_file;
  const unsigned int curr_time =
      std::chrono::system_clock::now().time_since_epoch().count();
  const std::string timing_stats_fname =
      base_name + "_timing_" + std::to_string(curr_time) + ".csv";
  timing_stats_file.open(timing_stats_fname, std::ios_base::app);
  timing_stats_file << "model,K,num_latent_factors,avg_full_dot_products,"
                       "preproc_time,comp_time"
                    << std::endl;
  const std::string timing_stats =
      (boost::format("%1%,%2%,%3%,%4%,%5%,%6%") % base_name % K %
       num_latent_factors % avg_full_dot_products % preprocess_time_s.count() %
       comp_time_s.count())
          .str();
  timing_stats_file << timing_stats << std::endl;
  timing_stats_file.close();
  return 0;
}
