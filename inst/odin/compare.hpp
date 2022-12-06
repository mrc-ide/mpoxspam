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
real_type ll_binom(real_type data_a, real_type data_b, real_type model_a, real_type model_b,
                   real_type exp_noise, rng_state_type& rng_state) {
  if (std::isnan(data_a) || std::isnan(data_b)) {
    return 0;
  }
  const real_type noise_a = dust::random::exponential<real_type>(rng_state, exp_noise);
  const real_type noise_b = dust::random::exponential<real_type>(rng_state, exp_noise);
  const real_type prob = exp_noise == std::numeric_limits<real_type>::infinity() && model_a == 0 && model_b == 0 ? 0 :
    (model_a + noise_a) / (model_a + model_b + noise_a + noise_b);
  return dust::density::binomial(data_a, data_a + data_b, prob, true);
}

// [[odin.dust::compare_data(Ytravel = real_type, Yendog = real_type, Yunk = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
compare(const typename T::real_type * state,
        const typename T::data_type& data,
        const typename T::internal_type internal,
        std::shared_ptr<const typename T::shared_type> shared,
        typename T::rng_state_type& rng_state) {
  typedef typename T::real_type real_type;
  const real_type newI = odin(newI);
  const real_type newIseed = odin(newIseed);
  const real_type time = odin(time);
  const real_type Yknown = data.Ytravel + data.Yendog;
  const real_type Y = Yknown + data.Yunk;

  real_type ret = 0;
  if (!std::isnan(Y)) {
    if (odin(use_nbinom)) {
      const real_type ll_cases = ll_nbinom(Y, newI, odin(kappa_cases), odin(exp_noise), rng_state);
      const real_type ll_travel = ll_binom(data.Ytravel, data.Yendog, newIseed, newI, odin(exp_noise), rng_state);
      ret = ll_cases + ll_travel;
    } else {
      const real_type delta = std::max(0.01, std::min(odin(delta1), odin(delta0) + odin(delta_slope) * time));
      constexpr real_type inf = std::numeric_limits<real_type>::infinity();

      const real_type t1 = newI < Y ? -inf : dust::density::binomial(Y, std::ceil(newI), delta, true);

      real_type t2 = 0;
      if (Yknown > 0) {
        const real_type newItotal = newIseed + newI;
        t2 = newItotal == 0 ? -inf : dust::density::binomial(std::ceil(data.Ytravel), std::ceil(Yknown), newIseed / newItotal, true);
      }
      ret = t1 + t2;
    }
  }
  return ret;
}
