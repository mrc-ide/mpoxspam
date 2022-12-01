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
  const real_type Y = dust::math::ceil(Yknown + data.Yunk);

  real_type ret = 0;
  if (!std::isnan(Y)) {
    const real_type delta = dust::math::max(static_cast<real_type>(0.01),
                                dust::math::min(odin(delta1), odin(delta0) + odin(delta_slope) * odin(time)));

    const real_type t1 = model_newI < Y ? -inf<real_type> : dust::density::binomial(Y, dust::math::ceil(model_newI), delta, true);

    real_type t2 = 0;
    if (Yknown > 0) {
      const real_type model_newItotal = model_newIseed + model_newI;
      t2 = model_newItotal == 0 ? -inf<real_type> : dust::density::binomial(dust::math::ceil(data.Ytravel), dust::math::ceil(Yknown), model_newIseed / model_newItotal, true);
    }
    ret = t1 + t2;
  }
  return ret;
}
