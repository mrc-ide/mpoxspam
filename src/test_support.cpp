#include <cpp11.hpp>
#include <dust/r/random.hpp>
#include <dust/random/cuda_compatibility.hpp>
#include "support.hpp"

[[cpp11::register]]
cpp11::sexp test_f(double x) {
  return cpp11::writable::doubles({f(x), fp(x), fpp(x), fppp(x)});
}

[[cpp11::register]]
cpp11::doubles test_g(double x) {
  return cpp11::writable::doubles({g(x), gp(x), gpp(x), gppp(x)});
}

[[cpp11::register]]
cpp11::doubles test_h(double x) {
  return cpp11::writable::doubles({h(x), hp(x), hpp(x), hppp(x)});
}

[[cpp11::register]]
double test_update_theta_vacc4_2(double theta_vacc, double amt_targetted) {
  return update_theta_vacc4_2(theta_vacc, amt_targetted);
}

[[cpp11::register]]
double test_ll_nbinom(double data, double model, double kappa,
                      double exp_noise, cpp11::sexp rng_ptr) {
  auto rng =
    dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  return ll_nbinom(data, model, kappa, exp_noise, rng->state(0));
}

[[cpp11::register]]
double test_ll_betabinom(double data_a, double data_b, double model_a, double model_b,
                         double rho, double exp_noise, cpp11::sexp rng_ptr) {
  auto rng =
    dust::random::r::rng_pointer_get<dust::random::xoshiro256plus>(rng_ptr);
  return ll_betabinom(data_a, data_b, model_a, model_b, rho, exp_noise, rng->state(0));
}
