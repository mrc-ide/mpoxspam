// There's a bit of a fight here with min, and none of this is going
// to be namespace-safe, so it's something that we should fix up
// later...
template <typename T>
__host__ __device__
T min(const T& a, const T&b) {
  return a < b ? a : b;
}

template <typename T>
__host__ __device__
T max(const T& a, const T&b) {
  return a > b ? a : b;
}

template <typename real_type>
constexpr real_type inf = std::numeric_limits<real_type>::infinity();

// [[odin.dust::compare_data(Ytravel = real_type)]]
// [[odin.dust::compare_data(Yendog = real_type)]]
// [[odin.dust::compare_data(Yunk = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
compare(const typename T::real_type * state,
        const typename T::data_type& data,
        const typename T::internal_type internal,
        std::shared_ptr<const typename T::shared_type> shared,
        typename T::rng_state_type& rng_state) {
  typedef typename T::real_type real_type;
  const real_type model_newI = odin(newI);
  const real_type model_newIseed = odin(newIseed);
  const real_type Yknown = data.Ytravel + data.Yendog;
  const real_type Y = std::ceil(Yknown + data.Yunk);

  real_type ret = 0;
  if (!std::isnan(Y)) {
    const real_type delta = max(static_cast<real_type>(0.01),
                                min(odin(delta1), odin(delta0) + odin(delta_slope) * odin(time)));

    const real_type t1 = model_newI < Y ? -inf<real_type> : dust::density::binomial(Y, std::ceil(model_newI), delta, true);

    real_type t2 = 0;
    if (Yknown > 0) {
      const real_type model_newItotal = model_newIseed + model_newI;
      t2 = model_newItotal == 0 ? -inf<real_type> : dust::density::binomial(std::ceil(data.Ytravel), std::ceil(Yknown), model_newIseed / model_newItotal, true);
    }
    ret = t1 + t2;
  }
  return ret;
}
