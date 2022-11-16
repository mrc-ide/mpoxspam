// Generated by dust (version 0.12.0) - do not edit
#include <cpp11.hpp>

[[cpp11::register]]
cpp11::sexp dust_m4_2_capabilities();

[[cpp11::register]]
cpp11::sexp dust_m4_2_gpu_info();
[[cpp11::register]]
SEXP dust_cpu_m4_2_alloc(cpp11::list r_pars, bool pars_multi, size_t time,
                         cpp11::sexp r_n_particles, size_t n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp gpu_config);

[[cpp11::register]]
SEXP dust_cpu_m4_2_run(SEXP ptr, size_t time_end);

[[cpp11::register]]
SEXP dust_cpu_m4_2_simulate(SEXP ptr, cpp11::sexp time_end);

[[cpp11::register]]
SEXP dust_cpu_m4_2_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_m4_2_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                SEXP r_time, SEXP r_set_initial_state);

[[cpp11::register]]
SEXP dust_cpu_m4_2_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
size_t dust_cpu_m4_2_time(SEXP ptr);

[[cpp11::register]]
void dust_cpu_m4_2_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_m4_2_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_cpu_m4_2_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_cpu_m4_2_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_cpu_m4_2_set_data(SEXP ptr, cpp11::list data, bool shared);

[[cpp11::register]]
SEXP dust_cpu_m4_2_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_m4_2_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood);

[[cpp11::register]]
void dust_cpu_m4_2_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_cpu_m4_2_n_state(SEXP ptr);
#include <dust/r/dust.hpp>

// Generated by odin.dust (version 0.2.26) - do not edit
template <typename real_type, typename T, typename U>
__host__ __device__ real_type fmodr(T x, U y) {
  real_type tmp = std::fmod(static_cast<real_type>(x),
                            static_cast<real_type>(y));
  if (tmp * y < 0) {
    tmp += y;
  }
  return tmp;
}

// These exist to support the model on the gpu, as in C++14 std::min
// and std::max are constexpr and error without --expt-relaxed-constexpr
template <typename T>
__host__ __device__ T odin_min(T x, T y) {
  return x < y ? x : y;
}

template <typename T>
__host__ __device__ T odin_max(T x, T y) {
  return x > y ? x : y;
}
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
// [[dust::class(m4_2)]]
// [[dust::param(N, has_default = TRUE, default_value = 750000L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta0, has_default = TRUE, default_value = 2.25, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta_freq, has_default = TRUE, default_value = 7L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta_sd, has_default = TRUE, default_value = 0.15, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(cumulative_partners_days, has_default = TRUE, default_value = 90L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta0, has_default = TRUE, default_value = 0.5, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta1, has_default = TRUE, default_value = 0.5, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_slope, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(etaf, has_default = TRUE, default_value = 0.005, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(etag, has_default = TRUE, default_value = 0.01, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma0, has_default = TRUE, default_value = 0.125, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma1, has_default = TRUE, default_value = 0.25, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(i0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(seedrate0, has_default = TRUE, default_value = 0.75, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(seedrate_sd, has_default = TRUE, default_value = 0.75, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_duration, has_default = TRUE, default_value = 55L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_freq, has_default = TRUE, default_value = 1L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_start_day, has_default = TRUE, default_value = 91L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_targetted, has_default = TRUE, default_value = 0.8, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class m4_2 {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    real_type N;
    real_type amt_random;
    real_type amt_targetted;
    real_type beta0;
    real_type beta_freq;
    real_type beta_sd;
    real_type cumulative_partners_days;
    real_type delta0;
    real_type delta1;
    real_type delta_slope;
    real_type etaf;
    real_type etag;
    real_type gamma0;
    real_type gamma1;
    real_type i0;
    real_type initial_E;
    real_type initial_Eseed;
    real_type initial_I;
    real_type initial_MEf;
    real_type initial_MEg;
    real_type initial_MEh;
    real_type initial_MIf;
    real_type initial_MIg;
    real_type initial_MIh;
    real_type initial_MSEf;
    real_type initial_MSEg;
    real_type initial_MSIf;
    real_type initial_MSIg;
    real_type initial_MSSf;
    real_type initial_MSSg;
    real_type initial_S_vacc;
    real_type initial_beta;
    real_type initial_cumulative_partners;
    real_type initial_cutf;
    real_type initial_cutg;
    real_type initial_cuth;
    real_type initial_cuts;
    real_type initial_dseedrate;
    real_type initial_newI;
    real_type initial_newIseed;
    real_type initial_seedrate;
    real_type initial_theta_vacc;
    real_type initial_thetaf;
    real_type initial_thetag;
    real_type initial_thetah;
    real_type seedrate0;
    real_type seedrate_sd;
    real_type vacc_amt;
    real_type vacc_duration;
    real_type vacc_fin_day;
    real_type vacc_freq;
    real_type vacc_start_day;
    real_type vacc_targetted;
    real_type xinit;
  };
  struct internal_type {
  };
  m4_2(const dust::pars_type<m4_2>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() {
    return 30;
  }
  std::vector<real_type> initial(size_t step) {
    std::vector<real_type> state(30);
    state[0] = shared->initial_thetaf;
    state[1] = shared->initial_MSEf;
    state[2] = shared->initial_MEf;
    state[3] = shared->initial_MSSf;
    state[4] = shared->initial_MSIf;
    state[5] = shared->initial_MIf;
    state[6] = shared->initial_thetag;
    state[7] = shared->initial_MSEg;
    state[8] = shared->initial_MEg;
    state[9] = shared->initial_MSSg;
    state[10] = shared->initial_MSIg;
    state[11] = shared->initial_MIg;
    state[12] = shared->initial_thetah;
    state[13] = shared->initial_MEh;
    state[14] = shared->initial_MIh;
    state[15] = shared->initial_E;
    state[16] = shared->initial_I;
    state[17] = shared->initial_newI;
    state[18] = shared->initial_Eseed;
    state[19] = shared->initial_newIseed;
    state[20] = shared->initial_cutf;
    state[21] = shared->initial_cutg;
    state[22] = shared->initial_cuth;
    state[23] = shared->initial_cuts;
    state[24] = shared->initial_seedrate;
    state[25] = shared->initial_dseedrate;
    state[26] = shared->initial_theta_vacc;
    state[27] = shared->initial_S_vacc;
    state[28] = shared->initial_beta;
    state[29] = shared->initial_cumulative_partners;
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type thetaf = state[0];
    const real_type MSEf = state[1];
    const real_type MEf = state[2];
    const real_type MSSf = state[3];
    const real_type MSIf = state[4];
    const real_type MIf = state[5];
    const real_type thetag = state[6];
    const real_type MSEg = state[7];
    const real_type MEg = state[8];
    const real_type MSSg = state[9];
    const real_type MSIg = state[10];
    const real_type MIg = state[11];
    const real_type thetah = state[12];
    const real_type MEh = state[13];
    const real_type MIh = state[14];
    const real_type E = state[15];
    const real_type I = state[16];
    const real_type newI = state[17];
    const real_type Eseed = state[18];
    const real_type newIseed = state[19];
    const real_type cutf = state[20];
    const real_type cutg = state[21];
    const real_type cuth = state[22];
    const real_type cuts = state[23];
    const real_type seedrate = state[24];
    const real_type dseedrate = state[25];
    const real_type theta_vacc = state[26];
    const real_type S_vacc = state[27];
    const real_type beta = state[28];
    real_type I_next = std::max(static_cast<real_type>(0), I + shared->gamma0 * E - shared->gamma1 * I);
    real_type dMIf = - shared->gamma1 * MIf + shared->gamma0 * MEf;
    real_type dMIg = - shared->gamma1 * MIg + shared->gamma0 * MEg;
    real_type dMIh = - shared->gamma1 * MIh + shared->gamma0 * MEh;
    real_type newI_next = newI + shared->gamma0 * E;
    real_type newIseed_next = newIseed + shared->gamma0 * Eseed;
    real_type time = step;
    state_next[26] = theta_vacc;
    real_type add_vaccine = (time >= shared->vacc_start_day) && (time <= shared->vacc_fin_day) && ((fmodr<real_type>((time - shared->vacc_start_day), shared->vacc_freq)) == 0);
    real_type beta_next = (fmodr<real_type>(time, shared->beta_freq) == 0 ? std::max(static_cast<real_type>(0), dust::random::normal<real_type>(rng_state, beta, shared->beta_sd)) : beta);
    real_type dseedrate_next = (fmodr<real_type>(time, shared->beta_freq) == 0 ? dust::random::normal<real_type>(rng_state, dseedrate, shared->seedrate_sd) : dseedrate);
    state_next[16] = I_next;
    state_next[5] = std::max(static_cast<real_type>(0), MIf + dMIf);
    state_next[11] = std::max(static_cast<real_type>(0), MIg + dMIg);
    state_next[14] = std::max(static_cast<real_type>(0), MIh + dMIh);
    state_next[17] = newI_next;
    state_next[19] = newIseed_next;
    real_type S_vacc_use = S_vacc * (1 - shared->amt_random);
    real_type rf = beta_next * 1.5 / (real_type) 7;
    real_type rg = beta_next * 1 / (real_type) 7;
    real_type seedrate_next = std::max(static_cast<real_type>(0), seedrate + dseedrate_next);
    real_type theta_vacc_use = update_theta_vacc4_2(theta_vacc, shared->amt_targetted);
    state_next[28] = beta_next;
    state_next[25] = dseedrate_next;
    real_type dot_thetaf = thetaf * theta_vacc_use;
    real_type dot_thetag = thetag * theta_vacc_use;
    real_type dot_thetah = thetah * theta_vacc_use;
    real_type red_f = (theta_vacc_use * fp(theta_vacc_use)) / (real_type) (theta_vacc * fp(theta_vacc));
    real_type red_g = (theta_vacc_use * gp(theta_vacc_use)) / (real_type) (theta_vacc * gp(theta_vacc));
    real_type transmseed = dust::random::poisson<real_type>(rng_state, seedrate_next);
    real_type trateh = std::max(static_cast<real_type>(0), beta_next * MIh * shared->N * hp(1) * S_vacc_use);
    state_next[27] = S_vacc_use;
    state_next[24] = seedrate_next;
    real_type Eseed_next = std::max(static_cast<real_type>(0), Eseed + transmseed - shared->gamma0 * Eseed);
    real_type MSf = dot_thetaf * S_vacc_use * fp(dot_thetaf) / (real_type) fp(1);
    real_type MSg = dot_thetag * S_vacc_use * gp(dot_thetag) / (real_type) gp(1);
    real_type meanfield_delta_si_f = (dot_thetaf * fpp(dot_thetaf) / (real_type) fp(dot_thetaf));
    real_type meanfield_delta_si_g = (dot_thetag * gpp(dot_thetag) / (real_type) gp(dot_thetag));
    real_type meanfield_delta_si_h = (1 + dot_thetah * hpp(dot_thetah) / (real_type) hp(dot_thetah));
    real_type transmh = dust::random::poisson<real_type>(rng_state, trateh);
    real_type u2f = (dot_thetaf * fpp(dot_thetaf) + std::pow(dot_thetaf, 2) * fppp(dot_thetaf)) / (real_type) fp(dot_thetaf);
    real_type u2g = (dot_thetag * gpp(dot_thetag) + std::pow(dot_thetag, 2) * gppp(dot_thetag)) / (real_type) gp(dot_thetag);
    state_next[23] = cuts + transmseed;
    real_type vaccine_scale_f = (add_vaccine ? (1 - shared->amt_random) * red_f : 1);
    real_type vaccine_scale_g = (add_vaccine ? (1 - shared->amt_random) * red_g : 1);
    real_type MSEf_vacc = MSEf * vaccine_scale_f;
    real_type MSEg_vacc = MSEg * vaccine_scale_g;
    real_type MSIf_vacc = MSIf * vaccine_scale_f;
    real_type MSIg_vacc = MSIg * vaccine_scale_g;
    real_type MSSf_vacc = MSSf * std::pow(vaccine_scale_f, 2);
    real_type MSSg_vacc = MSSg * std::pow(vaccine_scale_g, 2);
    real_type dthetah = - dot_thetah * (transmh + transmseed) / (real_type) (shared->N * hp(1));
    real_type u1f = meanfield_delta_si_f;
    real_type u1g = meanfield_delta_si_g;
    real_type u1h = meanfield_delta_si_h;
    state_next[18] = Eseed_next;
    state_next[22] = cuth + transmh;
    real_type dSh = hp(dot_thetah) * dthetah;
    real_type tratef = std::max(static_cast<real_type>(0), MSIf_vacc * shared->N * fp(1) * rf);
    real_type trateg = std::max(static_cast<real_type>(0), MSIg_vacc * shared->N * gp(1) * rg);
    real_type u2h = hppp(dot_thetah) * std::pow(dot_thetah, 2) / (real_type) hp(dot_thetah) + 2 * dot_thetah * hpp(dot_thetah) / (real_type) hp(dot_thetah) + u1h;
    state_next[12] = std::max(1.0000000000000001e-09, std::min(static_cast<real_type>(1), thetah + dthetah));
    real_type vf = u2f - std::pow(u1f, 2);
    real_type vg = u2g - std::pow(u1g, 2);
    real_type transmf = dust::random::poisson<real_type>(rng_state, tratef);
    real_type transmg = dust::random::poisson<real_type>(rng_state, trateg);
    real_type vh = u2h - std::pow(u1h, 2);
    real_type delta_si_f = ((transmf == 0 ? 0 : dust::random::normal<real_type>(rng_state, meanfield_delta_si_f, std::sqrt(vf / (real_type) transmf))));
    real_type delta_si_g = ((transmg == 0 ? 0 : dust::random::normal<real_type>(rng_state, meanfield_delta_si_g, std::sqrt(vg / (real_type) transmg))));
    real_type delta_si_h = ((transmh == 0 ? 0 : dust::random::normal<real_type>(rng_state, meanfield_delta_si_h, std::sqrt(vh / (real_type) transmh))));
    real_type dthetaf = - dot_thetaf * transmf / (real_type) (MSf * shared->N * fp(1));
    real_type dthetag = - dot_thetag * transmg / (real_type) (MSg * shared->N * gp(1));
    real_type newE = transmf + transmg + transmh + transmseed;
    real_type tauf = (transmf / (real_type) (transmf + transmg + transmh + transmseed));
    real_type taug = (transmg / (real_type) (transmf + transmg + transmh + transmseed));
    real_type tauh = ((transmh + transmseed) / (real_type) (transmf + transmg + transmh + transmseed));
    state_next[20] = cutf + transmf;
    state_next[21] = cutg + transmg;
    real_type E_next = std::max(static_cast<real_type>(0), E + newE - shared->gamma0 * E);
    real_type cumulative_partners_next = (1 + shared->etaf * shared->cumulative_partners_days) * (tauf * meanfield_delta_si_f + taug * dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) + tauh * dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf)) + (1 + shared->etag * shared->cumulative_partners_days) * (tauf * dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) + taug * meanfield_delta_si_g + tauh * dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag)) + shared->cumulative_partners_days * (tauf * dot_thetah * hp(dot_thetah) / (real_type) h(dot_thetah) + taug * dot_thetah * hp(dot_thetah) / (real_type) h(dot_thetah) + tauh * meanfield_delta_si_h);
    real_type dSf = fp(dot_thetaf) * dthetaf;
    real_type dSg = gp(dot_thetag) * dthetag;
    state_next[0] = std::max(1.0000000000000001e-09, std::min(static_cast<real_type>(1), thetaf + dthetaf));
    state_next[6] = std::max(1.0000000000000001e-09, std::min(static_cast<real_type>(1), thetag + dthetag));
    real_type dMEf = - shared->gamma0 * MEf + (- dSf) * (delta_si_f / (real_type) fp(1)) + ((- dSg) + (- dSh)) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1));
    real_type dMEg = - shared->gamma0 * MEg + (- dSg) * (delta_si_g / (real_type) gp(1)) + ((- dSf) + (- dSh)) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1));
    real_type dMEh = - shared->gamma0 * MEh + (- dSh) * (delta_si_h / (real_type) hp(1)) + ((- dSf) + (- dSg)) * (dot_thetah * hp(dot_thetah) / (real_type) h(dot_thetah) / (real_type) hp(1));
    real_type dMSEf = - shared->gamma1 * MSEf + 2 * shared->etaf * MSf * MEf - shared->etaf * MSEf + (- dSf) * (delta_si_f / (real_type) fp(1)) * (MSSf_vacc / (real_type) MSf) + (- dSg) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1)) * (MSSf_vacc / (real_type) MSf) + (- dSh) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1)) * (MSSf_vacc / (real_type) MSf);
    real_type dMSEg = - shared->gamma0 * MSEg + 2 * shared->etag * MSg * MEg - shared->etag * MSEg + (- dSg) * (delta_si_g / (real_type) gp(1)) * (MSSg_vacc / (real_type) MSg) + (- dSf) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1)) * (MSSg_vacc / (real_type) MSg) + (- dSh) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1)) * (MSSg_vacc / (real_type) MSg);
    real_type dMSIf = - rf * MSIf - shared->gamma1 * MSIf + shared->gamma0 * MSEf + 2 * shared->etaf * MSf * MIf - shared->etaf * MSIf + (- dSf) * (delta_si_f / (real_type) fp(1)) * (- MSIf_vacc / (real_type) MSf) + (- dSg) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1)) * (- MSIf_vacc / (real_type) MSf) + (- dSh) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1)) * (- MSIf_vacc / (real_type) MSf);
    real_type dMSIg = - rg * MSIg - shared->gamma1 * MSIg + shared->gamma0 * MSEg + 2 * shared->etag * MSg * MIg - shared->etag * MSIg + (- dSg) * (delta_si_g / (real_type) gp(1)) * (- MSIg_vacc / (real_type) MSg) + (- dSf) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1)) * (- MSIg_vacc / (real_type) MSg) + (- dSh) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1)) * (- MSIg_vacc / (real_type) MSg);
    real_type dMSSf = 1 * shared->etaf * std::pow(MSf, 2) - shared->etaf * MSSf_vacc - (- dSf) * (delta_si_f / (real_type) fp(1)) * MSSf_vacc / (real_type) MSf - ((- dSg) + (- dSh)) * (dot_thetaf * fp(dot_thetaf) / (real_type) f(dot_thetaf) / (real_type) fp(1)) * MSSf_vacc / (real_type) MSf;
    real_type dMSSg = 1 * shared->etag * std::pow(MSg, 2) - shared->etag * MSSg_vacc - (- dSg) * (delta_si_g / (real_type) gp(1)) * MSSg_vacc / (real_type) MSg - ((- dSf) + (- dSh)) * (dot_thetag * gp(dot_thetag) / (real_type) g(dot_thetag) / (real_type) gp(1)) * MSSg_vacc / (real_type) MSg;
    state_next[15] = E_next;
    state_next[29] = cumulative_partners_next;
    state_next[2] = std::max(static_cast<real_type>(0), MEf + dMEf);
    state_next[8] = std::max(static_cast<real_type>(0), MEg + dMEg);
    state_next[13] = std::max(static_cast<real_type>(0), MEh + dMEh);
    state_next[1] = std::max(static_cast<real_type>(0), MSEf_vacc + dMSEf);
    state_next[7] = std::max(static_cast<real_type>(0), MSEg_vacc + dMSEg);
    state_next[4] = std::max(static_cast<real_type>(0), MSIf_vacc + dMSIf);
    state_next[10] = std::max(static_cast<real_type>(0), MSIg_vacc + dMSIg);
    state_next[3] = std::max(static_cast<real_type>(0), MSSf_vacc + dMSSf);
    state_next[9] = std::max(static_cast<real_type>(0), MSSg_vacc + dMSSg);
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
__host__ __device__
real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace dust {
template<>
dust::pars_type<m4_2> dust_pars<m4_2>(cpp11::list user) {
  using real_type = typename m4_2::real_type;
  auto shared = std::make_shared<m4_2::shared_type>();
  m4_2::internal_type internal;
  shared->initial_Eseed = 0;
  shared->initial_MEf = 0;
  shared->initial_MEg = 0;
  shared->initial_MIf = 0;
  shared->initial_MIg = 0;
  shared->initial_MSEf = 0;
  shared->initial_MSEg = 0;
  shared->initial_MSIf = 0;
  shared->initial_MSIg = 0;
  shared->initial_MSSf = 1;
  shared->initial_MSSg = 1;
  shared->initial_S_vacc = 1;
  shared->initial_cumulative_partners = 0;
  shared->initial_cutf = 0;
  shared->initial_cutg = 0;
  shared->initial_cuth = 0;
  shared->initial_cuts = 0;
  shared->initial_dseedrate = 0;
  shared->initial_newI = 0;
  shared->initial_newIseed = 0;
  shared->initial_theta_vacc = 1;
  shared->initial_thetaf = 1;
  shared->initial_thetag = 1;
  shared->N = 750000;
  shared->beta0 = 2.25;
  shared->beta_freq = 7;
  shared->beta_sd = 0.14999999999999999;
  shared->cumulative_partners_days = 90;
  shared->delta0 = 0.5;
  shared->delta1 = 0.5;
  shared->delta_slope = 0;
  shared->etaf = 0.0050000000000000001;
  shared->etag = 0.01;
  shared->gamma0 = 0.125;
  shared->gamma1 = 0.25;
  shared->i0 = 0;
  shared->seedrate0 = 0.75;
  shared->seedrate_sd = 0.75;
  shared->vacc_duration = 55;
  shared->vacc_freq = 1;
  shared->vacc_start_day = 91;
  shared->vacc_targetted = 0.80000000000000004;
  shared->N = user_get_scalar<real_type>(user, "N", shared->N, NA_REAL, NA_REAL);
  shared->beta0 = user_get_scalar<real_type>(user, "beta0", shared->beta0, NA_REAL, NA_REAL);
  shared->beta_freq = user_get_scalar<real_type>(user, "beta_freq", shared->beta_freq, NA_REAL, NA_REAL);
  shared->beta_sd = user_get_scalar<real_type>(user, "beta_sd", shared->beta_sd, NA_REAL, NA_REAL);
  shared->cumulative_partners_days = user_get_scalar<real_type>(user, "cumulative_partners_days", shared->cumulative_partners_days, NA_REAL, NA_REAL);
  shared->delta0 = user_get_scalar<real_type>(user, "delta0", shared->delta0, NA_REAL, NA_REAL);
  shared->delta1 = user_get_scalar<real_type>(user, "delta1", shared->delta1, NA_REAL, NA_REAL);
  shared->delta_slope = user_get_scalar<real_type>(user, "delta_slope", shared->delta_slope, NA_REAL, NA_REAL);
  shared->etaf = user_get_scalar<real_type>(user, "etaf", shared->etaf, NA_REAL, NA_REAL);
  shared->etag = user_get_scalar<real_type>(user, "etag", shared->etag, NA_REAL, NA_REAL);
  shared->gamma0 = user_get_scalar<real_type>(user, "gamma0", shared->gamma0, NA_REAL, NA_REAL);
  shared->gamma1 = user_get_scalar<real_type>(user, "gamma1", shared->gamma1, NA_REAL, NA_REAL);
  shared->i0 = user_get_scalar<real_type>(user, "i0", shared->i0, NA_REAL, NA_REAL);
  shared->seedrate0 = user_get_scalar<real_type>(user, "seedrate0", shared->seedrate0, NA_REAL, NA_REAL);
  shared->seedrate_sd = user_get_scalar<real_type>(user, "seedrate_sd", shared->seedrate_sd, NA_REAL, NA_REAL);
  shared->vacc_duration = user_get_scalar<real_type>(user, "vacc_duration", shared->vacc_duration, NA_REAL, NA_REAL);
  shared->vacc_freq = user_get_scalar<real_type>(user, "vacc_freq", shared->vacc_freq, NA_REAL, NA_REAL);
  shared->vacc_start_day = user_get_scalar<real_type>(user, "vacc_start_day", shared->vacc_start_day, NA_REAL, NA_REAL);
  shared->vacc_targetted = user_get_scalar<real_type>(user, "vacc_targetted", shared->vacc_targetted, NA_REAL, NA_REAL);
  shared->initial_beta = shared->beta0;
  shared->initial_seedrate = shared->seedrate0;
  shared->vacc_amt = 0.65000000000000002 * 50000 / (real_type) shared->vacc_duration;
  shared->vacc_fin_day = shared->vacc_start_day + shared->vacc_duration;
  shared->xinit = shared->i0 / (real_type) shared->N;
  shared->amt_random = (1 - shared->vacc_targetted) * shared->vacc_amt / (real_type) shared->N;
  shared->amt_targetted = shared->vacc_targetted * shared->vacc_amt / (real_type) shared->N;
  shared->initial_E = shared->xinit * shared->N / (real_type) 2;
  shared->initial_I = shared->xinit * shared->N / (real_type) 2;
  shared->initial_MEh = shared->xinit / (real_type) 2;
  shared->initial_MIh = shared->xinit / (real_type) 2;
  shared->initial_thetah = 1 - shared->xinit;
  return dust::pars_type<m4_2>(shared, internal);
}
template <>
cpp11::sexp dust_info<m4_2>(const dust::pars_type<m4_2>& pars) {
  const m4_2::internal_type internal = pars.internal;
  const std::shared_ptr<const m4_2::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"thetaf", "MSEf", "MEf", "MSSf", "MSIf", "MIf", "thetag", "MSEg", "MEg", "MSSg", "MSIg", "MIg", "thetah", "MEh", "MIh", "E", "I", "newI", "Eseed", "newIseed", "cutf", "cutg", "cuth", "cuts", "seedrate", "dseedrate", "theta_vacc", "S_vacc", "beta", "cumulative_partners"});
  cpp11::writable::list dim(30);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({1});
  dim[6] = cpp11::writable::integers({1});
  dim[7] = cpp11::writable::integers({1});
  dim[8] = cpp11::writable::integers({1});
  dim[9] = cpp11::writable::integers({1});
  dim[10] = cpp11::writable::integers({1});
  dim[11] = cpp11::writable::integers({1});
  dim[12] = cpp11::writable::integers({1});
  dim[13] = cpp11::writable::integers({1});
  dim[14] = cpp11::writable::integers({1});
  dim[15] = cpp11::writable::integers({1});
  dim[16] = cpp11::writable::integers({1});
  dim[17] = cpp11::writable::integers({1});
  dim[18] = cpp11::writable::integers({1});
  dim[19] = cpp11::writable::integers({1});
  dim[20] = cpp11::writable::integers({1});
  dim[21] = cpp11::writable::integers({1});
  dim[22] = cpp11::writable::integers({1});
  dim[23] = cpp11::writable::integers({1});
  dim[24] = cpp11::writable::integers({1});
  dim[25] = cpp11::writable::integers({1});
  dim[26] = cpp11::writable::integers({1});
  dim[27] = cpp11::writable::integers({1});
  dim[28] = cpp11::writable::integers({1});
  dim[29] = cpp11::writable::integers({1});
  dim.names() = nms;
  cpp11::writable::list index(30);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = cpp11::writable::integers({6});
  index[6] = cpp11::writable::integers({7});
  index[7] = cpp11::writable::integers({8});
  index[8] = cpp11::writable::integers({9});
  index[9] = cpp11::writable::integers({10});
  index[10] = cpp11::writable::integers({11});
  index[11] = cpp11::writable::integers({12});
  index[12] = cpp11::writable::integers({13});
  index[13] = cpp11::writable::integers({14});
  index[14] = cpp11::writable::integers({15});
  index[15] = cpp11::writable::integers({16});
  index[16] = cpp11::writable::integers({17});
  index[17] = cpp11::writable::integers({18});
  index[18] = cpp11::writable::integers({19});
  index[19] = cpp11::writable::integers({20});
  index[20] = cpp11::writable::integers({21});
  index[21] = cpp11::writable::integers({22});
  index[22] = cpp11::writable::integers({23});
  index[23] = cpp11::writable::integers({24});
  index[24] = cpp11::writable::integers({25});
  index[25] = cpp11::writable::integers({26});
  index[26] = cpp11::writable::integers({27});
  index[27] = cpp11::writable::integers({28});
  index[28] = cpp11::writable::integers({29});
  index[29] = cpp11::writable::integers({30});
  index.names() = nms;
  size_t len = 30;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

cpp11::sexp dust_m4_2_capabilities() {
  return dust::r::dust_capabilities<m4_2>();
}

cpp11::sexp dust_m4_2_gpu_info() {
  return dust::gpu::r::gpu_info();
}
using model_cpu = dust::dust_cpu<m4_2>;

SEXP dust_cpu_m4_2_alloc(cpp11::list r_pars, bool pars_multi, size_t time,
                             cpp11::sexp r_n_particles, size_t n_threads,
                             cpp11::sexp r_seed, bool deterministic,
                             cpp11::sexp gpu_config) {
  return dust::r::dust_cpu_alloc<m4_2>(r_pars, pars_multi, time, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        gpu_config);
}

SEXP dust_cpu_m4_2_run(SEXP ptr, size_t time_end) {
  return dust::r::dust_run<model_cpu>(ptr, time_end);
}

SEXP dust_cpu_m4_2_simulate(SEXP ptr, cpp11::sexp time_end) {
  return dust::r::dust_simulate<model_cpu>(ptr, time_end);
}

SEXP dust_cpu_m4_2_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model_cpu>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_cpu_m4_2_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                SEXP r_time, SEXP r_set_initial_state) {
  return dust::r::dust_update_state<model_cpu>(ptr, r_pars, r_state, r_time,
                                               r_set_initial_state);
}

SEXP dust_cpu_m4_2_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model_cpu>(ptr, r_index);
}

size_t dust_cpu_m4_2_time(SEXP ptr) {
  return dust::r::dust_time<model_cpu>(ptr);
}

void dust_cpu_m4_2_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_m4_2_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model_cpu>(ptr, r_weights);
}

SEXP dust_cpu_m4_2_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model_cpu>(ptr, first_only, last_only);
}

SEXP dust_cpu_m4_2_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model_cpu>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_cpu_m4_2_set_data(SEXP ptr, cpp11::list data,
                                       bool shared) {
  dust::r::dust_set_data<model_cpu>(ptr, data, shared);
  return R_NilValue;
}

SEXP dust_cpu_m4_2_compare_data(SEXP ptr) {
  return dust::r::dust_compare_data<model_cpu>(ptr);
}

SEXP dust_cpu_m4_2_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood) {
  return dust::r::dust_filter<model_cpu>(ptr, time_end,
                                                save_trajectories,
                                                time_snapshot,
                                                min_log_likelihood);
}

void dust_cpu_m4_2_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model_cpu>(ptr, n_threads);
}

int dust_cpu_m4_2_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model_cpu>(ptr);
}
