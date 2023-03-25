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
cpp11::doubles test_h(double x, double hs, double hr) {
  return cpp11::writable::doubles({h(x, hs, hr), hp(x, hs, hr), hpp(x, hs, hr), hppp(x, hs, hr)});
}

[[cpp11::register]]
cpp11::doubles test_hu(double x, double vr, double V1, double V2, double v1eff, double v2eff, double thetav, double hs, double hr) {
  return cpp11::writable::doubles({hu(x, vr, V1, V2, v1eff, v2eff, thetav, hs, hr),
                                   hup(x, vr, V1, V2, v1eff, v2eff, thetav, hs, hr),
                                   hupp(x, vr, V1, V2, v1eff, v2eff, thetav, hs, hr),
                                   huppp(x, vr, V1, V2, v1eff, v2eff, thetav, hs, hr)});
}

[[cpp11::register]]
double test_update_theta_vacc4_2(double theta_vacc, double amt_targetted, double hshape, double hrate) {
  return update_theta_vacc4_2(theta_vacc, amt_targetted, hshape, hrate);
}

[[cpp11::register]]
double test_update_theta_vacc4_3(double prop_vacc_targetted, double hshape, double hrate) {
  return update_theta_vacc4_3(prop_vacc_targetted, hshape, hrate);
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
