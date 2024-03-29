// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// model.cpp
cpp11::sexp dust_model_gpu_info();
extern "C" SEXP _mpoxspam_dust_model_gpu_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_model_gpu_info());
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time, cpp11::sexp r_n_particles, int n_threads, cpp11::sexp r_seed, bool deterministic, cpp11::sexp gpu_config, cpp11::sexp ode_control);
extern "C" SEXP _mpoxspam_dust_cpu_model_alloc(SEXP r_pars, SEXP pars_multi, SEXP r_time, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP deterministic, SEXP gpu_config, SEXP ode_control) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<bool>>(deterministic), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(gpu_config), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ode_control)));
  END_CPP11
}
// model.cpp
cpp11::sexp dust_cpu_model_capabilities();
extern "C" SEXP _mpoxspam_dust_cpu_model_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_capabilities());
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_run(SEXP ptr, cpp11::sexp r_time_end);
extern "C" SEXP _mpoxspam_dust_cpu_model_run(SEXP ptr, SEXP r_time_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time_end)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_simulate(SEXP ptr, cpp11::sexp time_end);
extern "C" SEXP _mpoxspam_dust_cpu_model_simulate(SEXP ptr, SEXP time_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(time_end)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _mpoxspam_dust_cpu_model_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size);
extern "C" SEXP _mpoxspam_dust_cpu_model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_update_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_pars), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_time), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_set_initial_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(index), cpp11::as_cpp<cpp11::decay_t<SEXP>>(reset_step_size)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _mpoxspam_dust_cpu_model_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_time(SEXP ptr);
extern "C" SEXP _mpoxspam_dust_cpu_model_time(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_time(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// model.cpp
void dust_cpu_model_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _mpoxspam_dust_cpu_model_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_cpu_model_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _mpoxspam_dust_cpu_model_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_rng_state(SEXP ptr, bool first_only, bool last_only);
extern "C" SEXP _mpoxspam_dust_cpu_model_rng_state(SEXP ptr, SEXP first_only, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(first_only), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _mpoxspam_dust_cpu_model_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_set_data(SEXP ptr, cpp11::list data, bool shared);
extern "C" SEXP _mpoxspam_dust_cpu_model_set_data(SEXP ptr, SEXP data, SEXP shared) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data), cpp11::as_cpp<cpp11::decay_t<bool>>(shared)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_compare_data(SEXP ptr);
extern "C" SEXP _mpoxspam_dust_cpu_model_compare_data(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_filter(SEXP ptr, SEXP time_end, bool save_trajectories, cpp11::sexp time_snapshot, cpp11::sexp min_log_likelihood);
extern "C" SEXP _mpoxspam_dust_cpu_model_filter(SEXP ptr, SEXP time_end, SEXP save_trajectories, SEXP time_snapshot, SEXP min_log_likelihood) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(time_end), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(time_snapshot), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(min_log_likelihood)));
  END_CPP11
}
// model.cpp
void dust_cpu_model_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _mpoxspam_dust_cpu_model_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_cpu_model_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// model.cpp
int dust_cpu_model_n_state(SEXP ptr);
extern "C" SEXP _mpoxspam_dust_cpu_model_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// model.cpp
void dust_cpu_model_set_stochastic_schedule(SEXP ptr, SEXP time);
extern "C" SEXP _mpoxspam_dust_cpu_model_set_stochastic_schedule(SEXP ptr, SEXP time) {
  BEGIN_CPP11
    dust_cpu_model_set_stochastic_schedule(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(time));
    return R_NilValue;
  END_CPP11
}
// model.cpp
SEXP dust_cpu_model_ode_statistics(SEXP ptr);
extern "C" SEXP _mpoxspam_dust_cpu_model_ode_statistics(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_model_ode_statistics(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// test_support.cpp
cpp11::sexp test_f(double x);
extern "C" SEXP _mpoxspam_test_f(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_f(cpp11::as_cpp<cpp11::decay_t<double>>(x)));
  END_CPP11
}
// test_support.cpp
cpp11::doubles test_g(double x);
extern "C" SEXP _mpoxspam_test_g(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_g(cpp11::as_cpp<cpp11::decay_t<double>>(x)));
  END_CPP11
}
// test_support.cpp
cpp11::doubles test_h(double x, double hs, double hr);
extern "C" SEXP _mpoxspam_test_h(SEXP x, SEXP hs, SEXP hr) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_h(cpp11::as_cpp<cpp11::decay_t<double>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(hs), cpp11::as_cpp<cpp11::decay_t<double>>(hr)));
  END_CPP11
}
// test_support.cpp
cpp11::doubles test_hu(double x, double vr, double V1, double V2, double v1eff, double v2eff, double thetav, double hs, double hr);
extern "C" SEXP _mpoxspam_test_hu(SEXP x, SEXP vr, SEXP V1, SEXP V2, SEXP v1eff, SEXP v2eff, SEXP thetav, SEXP hs, SEXP hr) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_hu(cpp11::as_cpp<cpp11::decay_t<double>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(vr), cpp11::as_cpp<cpp11::decay_t<double>>(V1), cpp11::as_cpp<cpp11::decay_t<double>>(V2), cpp11::as_cpp<cpp11::decay_t<double>>(v1eff), cpp11::as_cpp<cpp11::decay_t<double>>(v2eff), cpp11::as_cpp<cpp11::decay_t<double>>(thetav), cpp11::as_cpp<cpp11::decay_t<double>>(hs), cpp11::as_cpp<cpp11::decay_t<double>>(hr)));
  END_CPP11
}
// test_support.cpp
double test_update_theta_vacc4_2(double theta_vacc, double amt_targetted, double hshape, double hrate);
extern "C" SEXP _mpoxspam_test_update_theta_vacc4_2(SEXP theta_vacc, SEXP amt_targetted, SEXP hshape, SEXP hrate) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_update_theta_vacc4_2(cpp11::as_cpp<cpp11::decay_t<double>>(theta_vacc), cpp11::as_cpp<cpp11::decay_t<double>>(amt_targetted), cpp11::as_cpp<cpp11::decay_t<double>>(hshape), cpp11::as_cpp<cpp11::decay_t<double>>(hrate)));
  END_CPP11
}
// test_support.cpp
double test_update_theta_vacc4_3(double prop_vacc_targetted, double hshape, double hrate);
extern "C" SEXP _mpoxspam_test_update_theta_vacc4_3(SEXP prop_vacc_targetted, SEXP hshape, SEXP hrate) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_update_theta_vacc4_3(cpp11::as_cpp<cpp11::decay_t<double>>(prop_vacc_targetted), cpp11::as_cpp<cpp11::decay_t<double>>(hshape), cpp11::as_cpp<cpp11::decay_t<double>>(hrate)));
  END_CPP11
}
// test_support.cpp
double test_ll_nbinom(double data, double model, double kappa, double exp_noise, cpp11::sexp rng_ptr);
extern "C" SEXP _mpoxspam_test_ll_nbinom(SEXP data, SEXP model, SEXP kappa, SEXP exp_noise, SEXP rng_ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_ll_nbinom(cpp11::as_cpp<cpp11::decay_t<double>>(data), cpp11::as_cpp<cpp11::decay_t<double>>(model), cpp11::as_cpp<cpp11::decay_t<double>>(kappa), cpp11::as_cpp<cpp11::decay_t<double>>(exp_noise), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(rng_ptr)));
  END_CPP11
}
// test_support.cpp
double test_ll_betabinom(double data_a, double data_b, double model_a, double model_b, double rho, double exp_noise, cpp11::sexp rng_ptr);
extern "C" SEXP _mpoxspam_test_ll_betabinom(SEXP data_a, SEXP data_b, SEXP model_a, SEXP model_b, SEXP rho, SEXP exp_noise, SEXP rng_ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(test_ll_betabinom(cpp11::as_cpp<cpp11::decay_t<double>>(data_a), cpp11::as_cpp<cpp11::decay_t<double>>(data_b), cpp11::as_cpp<cpp11::decay_t<double>>(model_a), cpp11::as_cpp<cpp11::decay_t<double>>(model_b), cpp11::as_cpp<cpp11::decay_t<double>>(rho), cpp11::as_cpp<cpp11::decay_t<double>>(exp_noise), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(rng_ptr)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_mpoxspam_dust_cpu_model_alloc",                   (DL_FUNC) &_mpoxspam_dust_cpu_model_alloc,                   9},
    {"_mpoxspam_dust_cpu_model_capabilities",            (DL_FUNC) &_mpoxspam_dust_cpu_model_capabilities,            0},
    {"_mpoxspam_dust_cpu_model_compare_data",            (DL_FUNC) &_mpoxspam_dust_cpu_model_compare_data,            1},
    {"_mpoxspam_dust_cpu_model_filter",                  (DL_FUNC) &_mpoxspam_dust_cpu_model_filter,                  5},
    {"_mpoxspam_dust_cpu_model_n_state",                 (DL_FUNC) &_mpoxspam_dust_cpu_model_n_state,                 1},
    {"_mpoxspam_dust_cpu_model_ode_statistics",          (DL_FUNC) &_mpoxspam_dust_cpu_model_ode_statistics,          1},
    {"_mpoxspam_dust_cpu_model_reorder",                 (DL_FUNC) &_mpoxspam_dust_cpu_model_reorder,                 2},
    {"_mpoxspam_dust_cpu_model_resample",                (DL_FUNC) &_mpoxspam_dust_cpu_model_resample,                2},
    {"_mpoxspam_dust_cpu_model_rng_state",               (DL_FUNC) &_mpoxspam_dust_cpu_model_rng_state,               3},
    {"_mpoxspam_dust_cpu_model_run",                     (DL_FUNC) &_mpoxspam_dust_cpu_model_run,                     2},
    {"_mpoxspam_dust_cpu_model_set_data",                (DL_FUNC) &_mpoxspam_dust_cpu_model_set_data,                3},
    {"_mpoxspam_dust_cpu_model_set_index",               (DL_FUNC) &_mpoxspam_dust_cpu_model_set_index,               2},
    {"_mpoxspam_dust_cpu_model_set_n_threads",           (DL_FUNC) &_mpoxspam_dust_cpu_model_set_n_threads,           2},
    {"_mpoxspam_dust_cpu_model_set_rng_state",           (DL_FUNC) &_mpoxspam_dust_cpu_model_set_rng_state,           2},
    {"_mpoxspam_dust_cpu_model_set_stochastic_schedule", (DL_FUNC) &_mpoxspam_dust_cpu_model_set_stochastic_schedule, 2},
    {"_mpoxspam_dust_cpu_model_simulate",                (DL_FUNC) &_mpoxspam_dust_cpu_model_simulate,                2},
    {"_mpoxspam_dust_cpu_model_state",                   (DL_FUNC) &_mpoxspam_dust_cpu_model_state,                   2},
    {"_mpoxspam_dust_cpu_model_time",                    (DL_FUNC) &_mpoxspam_dust_cpu_model_time,                    1},
    {"_mpoxspam_dust_cpu_model_update_state",            (DL_FUNC) &_mpoxspam_dust_cpu_model_update_state,            7},
    {"_mpoxspam_dust_model_gpu_info",                    (DL_FUNC) &_mpoxspam_dust_model_gpu_info,                    0},
    {"_mpoxspam_test_f",                                 (DL_FUNC) &_mpoxspam_test_f,                                 1},
    {"_mpoxspam_test_g",                                 (DL_FUNC) &_mpoxspam_test_g,                                 1},
    {"_mpoxspam_test_h",                                 (DL_FUNC) &_mpoxspam_test_h,                                 3},
    {"_mpoxspam_test_hu",                                (DL_FUNC) &_mpoxspam_test_hu,                                9},
    {"_mpoxspam_test_ll_betabinom",                      (DL_FUNC) &_mpoxspam_test_ll_betabinom,                      7},
    {"_mpoxspam_test_ll_nbinom",                         (DL_FUNC) &_mpoxspam_test_ll_nbinom,                         5},
    {"_mpoxspam_test_update_theta_vacc4_2",              (DL_FUNC) &_mpoxspam_test_update_theta_vacc4_2,              4},
    {"_mpoxspam_test_update_theta_vacc4_3",              (DL_FUNC) &_mpoxspam_test_update_theta_vacc4_3,              3},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_mpoxspam(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
