#include <cpp11.hpp>
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
