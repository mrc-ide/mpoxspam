#include "uniroot.hpp"

// [[odin.dust::register]]
inline double f(double x) {
  return (2809 + 2000 * x + 95 * x * x) / 4904.0;
}

// [[odin.dust::register]]
inline double fp(double x) {
  return (2000 + 2 * 95 * x) / 4904.0;
}

// [[odin.dust::register]]
inline double fpp(double x) {
  return (2.0 * 95.0) / 4904.0;
}

// [[odin.dust::register]]
inline double fppp(double x) {
  return 0;
}

// [[odin.dust::register]]
inline double g(double x) {
  return (2943 + 1009 * x + 477 * x * x + 475 * x * x * x) / 4904.0;
}

// [[odin.dust::register]]
inline double gp(double x) {
  return (1009 + 2 * 477 * x + 3 * 475 * x * x) / 4904.0;
}

// [[odin.dust::register]]
inline double gpp(double x) {
  return (2 * 477 + 2 * 3 * 475 * x) / 4904.0;
}

// [[odin.dust::register]]
inline double gppp(double x) {
  return (2 * 3 * 475.0) / 4904.0;
}

constexpr double hshape = 0.26;
constexpr double hrate = 1.85 * 7;

// [[odin.dust::register]]
inline double h(double x) {
  return std::pow(1 - std::log(x) / hrate, -hshape);
}

// [[odin.dust::register]]
inline double hp(double x) {
  return hshape * std::pow(1 - std::log(x) / hrate, -hshape - 1) / (hrate * x);
}

// [[odin.dust::register]]
inline double hpp(double x) {
  return hshape * std::pow((hrate - std::log(x)) / hrate, -hshape) *
    (-hrate + hshape + std::log(x) + 1) / (x * x * std::pow(hrate - std::log(x), 2));
}

// [[odin.dust::register]]
inline double hppp(double x) {
  return hshape * std::pow((hrate - std::log(x)) / hrate, -hshape) *
    (hshape * hshape + 3 * hshape + 2 * std::pow(hrate - std::log(x), 2) - 3 * (hrate - std::log(x)) * (hshape + 1) + 2) / (x * x * x * std::pow(hrate - std::log(x), 3));
}

// [[odin.dust::register]]
inline double update_theta_vacc4_2(double theta_vacc, double amt_targetted) {
  const double p0 = f(theta_vacc) * g(theta_vacc) * h(theta_vacc);
  const double p1 = std::max(0.01, p0 - amt_targetted);
  const auto fn = [&](double x) { return f(exp(x)) * g(exp(x)) * h(exp(x)) - p1; };
  return std::exp(uniroot_brent<double>(fn, -1e4, 0, 1e-6, 1000));
}
