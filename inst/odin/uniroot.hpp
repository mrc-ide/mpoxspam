#pragma once
#include <limits>
#include <stdexcept>

// This can be generic
template <typename real_type>
struct uniroot_result {
  real_type x;
  real_type fx;
  int iterations;
  bool converged;
};

// From zeroin.c, in brent.shar
template <typename real_type, typename F>
uniroot_result<real_type> uniroot_brent_run(F f, real_type a, real_type b,
                                            real_type tol, int max_iterations) {
  real_type fa = f(a);
  real_type fb = f(b);
  int iterations = 0;

  if (fa == 0) {
    return uniroot_result<real_type>{a, fa, iterations, true};
  }
  if (fb == 0) {
    return uniroot_result<real_type>{b, fb, iterations, true};
  }
  if (fa * fb > 0) {
    // Same sign
    //
    // Be careful with result, as it will contain junk
    return uniroot_result<real_type>{a, fa, iterations, false};
  }

  real_type c = a;
  real_type fc = fa;   // c = a, f(c) = f(a)
  const real_type eps = std::numeric_limits<real_type>::epsilon();


  for (; iterations < max_iterations; ++iterations) { // Main iteration loop
    // Distance from the last but one to the last approximation
    const real_type prev_step = b - a;

    // Interpolation step is calculated in the form p/q; division
    // operations is dlayed until the last moment
    real_type p;
    real_type q;

    if (std::abs(fc) < std::abs(fb)) {
      // Swap data for b to be the best approximation
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    // Actual tolerance
    const real_type tol_act = 2 * eps * std::abs(b) + tol / 2;
    // Step at this iteration
    real_type new_step = (c - b) / 2;

    if (std::abs(new_step) <= tol_act || fb == 0) {
      // Acceptable approximation is found
      return uniroot_result<real_type>{b, fb, iterations, true};
    }

    // Decide if the interpolation can be tried
    //
    // If prev_step was large enough and was in true direction, then
    // interpolation can be tried
    if (std::abs(prev_step) >= tol_act && std::abs(fa) > std::abs(fb)) {
      // interpolation
      const real_type cb = c - b;
      if (a == c) {
        // If we have only two distinct points linear interpolation
        // can only be applied
        const real_type t1 = fb / fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        // Quadric inverse interpolation
        q = fa / fc;
        const real_type t1 = fb / fc;
        const real_type t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > 0) {
	// p was calculated with the opposite sign; make p positive
	// and assign possible minus to q
        q = -q;
      } else {
        p = -p;
      }

      // If b + p / q falls in [b, c] and isn't too large it is
      // accepted
      //
      // If p / q is too large then the bissection procedure can
      // reduce [b,c] range to more extent
      if (p < (0.75 * cb * q - std::abs(tol_act * q) / 2) &&
          p < std::abs(prev_step * q / 2)) {
        new_step = p / q;
      }
    }

    // Adjust the step to be not less than tolerance
    if (std::abs(new_step) < tol_act) {
      if (new_step > 0) {
        new_step = tol_act;
      } else {
        new_step = -tol_act;
      }
    }

    // Save the previous approximation
    a = b;
    fa = fb;
    // Do step to a new approximation
    b += new_step;
    fb = f(b);
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      // Adjust c for it to have a sign opposite to that of b
      c = a;  fc = fa;
    }
  }

  return uniroot_result<real_type>{b, fb, iterations, false};
}

template <typename real_type, typename F>
real_type uniroot_brent(F f, real_type a, real_type b,
                        real_type tol, int max_iterations) {
  const auto result = uniroot_brent_run(f, a, b, tol, max_iterations);
  if (!result.converged) {
    throw std::runtime_error("some error");
  }
  return result.x;
}
