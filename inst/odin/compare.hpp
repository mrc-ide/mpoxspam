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
  const real_type Y = std::ceil(Yknown + data.Yunk);

  real_type ret = 0;
  if (!std::isnan(Y)) {
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
  return ret;
}
