// [[odin.dust::cpp_std(C++14)]]
// [[odin.dust::linking_to(lostturnip)]]
#include <lostturnip.hpp>
#include <dust/random/math.hpp>
#include <dust/random/density.hpp>
#include <dust/random/random.hpp>

template <typename real_type>
constexpr real_type inf = std::numeric_limits<real_type>::infinity();

template <typename real_type>
constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type f(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2809 + 2000 * x + 95 * x * x) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type fp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2000 + 2 * 95 * x) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type fpp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2 * 95) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type fppp(real_type x) {
  return 0;
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type g(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2943 + 1009 * x + 477 * x * x + 475 * x * x * x) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type gp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (1009 + 2 * 477 * x + 3 * 475 * x * x) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type gpp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2 * 477 + 2 * 3 * 475 * x) / static_cast<real_type>(4904);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type gppp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  return (2 * 3 * 475) / static_cast<real_type>(4904);
}

// There's a bit of a fight here with hrate and hshape because we need
// these to be in the correct precision given the their functions.
template <typename real_type>
constexpr real_type hshape = 0.26;
template <typename real_type>
constexpr real_type hrate = 1.85 * 7;

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type h(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  const real_type hs = hshape<real_type>;
  const real_type hr = hrate<real_type>;
  return dust::math::pow(1 - dust::math::log(x) / hr, -hs);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type hp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  const real_type hs = hshape<real_type>;
  const real_type hr = hrate<real_type>;
  return hs * dust::math::pow(1 - dust::math::log(x) / hr, -hs - 1) / (hr * x);
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type hpp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  const real_type hs = hshape<real_type>;
  const real_type hr = hrate<real_type>;
  return hs * dust::math::pow((hr - dust::math::log(x)) / hr, -hs) *
    (-hr + hs + dust::math::log(x) + 1) / (x * x * dust::math::pow(hr - dust::math::log(x), 2));
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type hppp(real_type x) {
  static_assert(std::is_floating_point<real_type>::value, "use with integral type");
  const real_type hs = hshape<real_type>;
  const real_type hr = hrate<real_type>;
  return hs * dust::math::pow((hr - dust::math::log(x)) / hr, -hs) *
    (hs * hs + 3 * hs + 2 * dust::math::pow(hr - dust::math::log(x), 2) - 3 * (hr - dust::math::log(x)) * (hs + 1) + 2) / (x * x * x * dust::math::pow(hr - dust::math::log(x), 3));
}

// [[odin.dust::register]]
template <typename real_type>
__host__ __device__
real_type update_theta_vacc4_2(real_type theta_vacc, real_type amt_targetted) {
  // Somewhat annoyingly, we can do this as constexpr in gcc but it's
  // not portable as sqrt() is not constexpr in any version of the
  // standard. It's possible that the compiler will work this out for
  // us?
  const real_type tol = dust::math::sqrt(eps<real_type>);
  const real_type p0 = f(theta_vacc) * g(theta_vacc) * h(theta_vacc);
  const real_type p1 = dust::math::max(static_cast<real_type>(0.01), p0 - amt_targetted);
  // There's a small optimisation that can be made here by avoiding
  // doing log(exp(x)) in h
  const auto fn = [&](real_type x) {
                    const auto exp_x = dust::math::exp(x);
                    return f(exp_x) * g(exp_x) * dust::math::pow(1 - x / hrate<real_type>, -hshape<real_type>) - p1;
                  };
  return dust::math::exp(lostturnip::find<real_type>(fn, -1e4, 0, tol, 1000));
}

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type ll_nbinom(real_type data, real_type model, real_type kappa,
                    real_type exp_noise, rng_state_type& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_type mu = model +
    dust::random::exponential<real_type>(rng_state, exp_noise);
  return dust::density::negative_binomial_mu(data, kappa, mu, true);
}

template <typename real_type, typename rng_state_type>
__host__ __device__
real_type ll_betabinom(real_type data_a, real_type data_b,
                       real_type model_a, real_type model_b,
                       real_type rho, real_type exp_noise,
                       rng_state_type& rng_state) {
  if (std::isnan(data_a) || std::isnan(data_b)) {
    return 0;
  }
  const real_type noise_a = dust::random::exponential<real_type>(rng_state, exp_noise);
  const real_type noise_b = dust::random::exponential<real_type>(rng_state, exp_noise);
  real_type prob_a = (model_a + noise_a) /
    (model_a + noise_a + model_b + noise_b);
  return dust::density::beta_binomial(data_a, data_a + data_b, prob_a, rho, true);
}

template <typename real_type>
__host__ __device__
real_type calc_delta(real_type delta0, real_type delta1, real_type delta_slope,
                    real_type time) {
  real_type delta_curr = delta0 + delta_slope * time;
  return dust::math::max(static_cast<real_type>(0.01),
                         dust::math::min(delta1, delta_curr));
}
