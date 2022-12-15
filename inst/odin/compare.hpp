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
  const real_type delta = calc_delta(odin(delta0), odin(delta1),
                                     odin(delta_slope), odin(time));
  const real_type cases = dust::math::ceil(model_newI  * delta);
  const real_type model_newIseed = odin(newIseed);
  const real_type Yknown = data.Ytravel + data.Yendog;
  const real_type Y = dust::math::ceil(Yknown + data.Yunk);

  real_type ret = 0;
  if (!std::isnan(Y)) {
    real_type ll_cases = ll_nbinom(Y, cases, odin(kappa_cases),
                                   odin(exp_noise), rng_state);
    real_type ll_travel = ll_betabinom(data.Ytravel, data.Yendog,
                                       model_newIseed, model_newI,
                                       odin(rho_travel), odin(exp_noise),
                                       rng_state);
    ret = ll_cases + ll_travel;
  }

  return ret;
}
