#include <gmpxx.h>

#include <complex>
#include <iostream>
#include <optional>
#include <vector>

const int PREC = 32;

// function prototype for fortran interface
extern "C" {
void thiele_pade_mp_api(int n_par, const std::complex<double> *x_ref,
                        const std::complex<double> *y_ref,
                        const std::complex<double> *x_query,
                        std::complex<double> *y_query, int num_query);
}

/// @brief  datatype of a complex number using arbitrary precision numbers for
///         real and imaginary part
class ComplexGMP {
public:
  /// real part
  mpf_class real_mp;
  /// imag part
  mpf_class imag_mp;

  /// @brief initialize a complex number 0 + 0i
  ComplexGMP() : real_mp(0), imag_mp(0) {}

  /// @brief initialize a complex number from real and imag part using arbitrary
  ///        precision numbers
  /// @param real_mp real part arbitrary precision
  /// @param imag_mp imaginary part arbitrary precision
  ComplexGMP(const mpf_class &real_mp, const mpf_class &imag_mp)
      : real_mp(real_mp), imag_mp(imag_mp) {}

  /// @brief initialize a complex number from a complex double
  /// @param z complex number
  ComplexGMP(std::complex<double> z) : real_mp(real(z)), imag_mp(imag(z)) {}

  /// @brief addition of two arbitrary precision complex numbers
  /// @param rhs right summand
  /// @return arbitrary precision comoplex number
  ComplexGMP operator+(const ComplexGMP &rhs) const {
    return ComplexGMP(real_mp + rhs.real_mp, imag_mp + rhs.imag_mp);
  }

  /// @brief substraction of two arbitrary precision complex numbers
  /// @param rhs subtrahend
  /// @return arbitrary precision comoplex number
  ComplexGMP operator-(const ComplexGMP &rhs) const {
    return ComplexGMP(real_mp - rhs.real_mp, imag_mp - rhs.imag_mp);
  }

  /// @brief multiplication of two arbitrary precision complex numbers
  /// @param rhs right factor
  /// @return arbitrary precision comoplex number
  ComplexGMP operator*(const ComplexGMP &rhs) const {
    mpf_class newreal_mp = real_mp * rhs.real_mp - imag_mp * rhs.imag_mp;
    mpf_class newimag_mp = real_mp * rhs.imag_mp + imag_mp * rhs.real_mp;
    return ComplexGMP(newreal_mp, newimag_mp);
  }

  /// @brief division of two arbitrary precision complex numbers
  /// @param rhs denominator
  /// @return arbitrary precision comoplex number
  ComplexGMP operator/(const ComplexGMP &rhs) const {
    mpf_class denominator =
        rhs.real_mp * rhs.real_mp + rhs.imag_mp * rhs.imag_mp;
    mpf_class newreal_mp =
        (real_mp * rhs.real_mp + rhs.imag_mp * imag_mp) / denominator;
    mpf_class newimag_mp =
        (rhs.real_mp * imag_mp - real_mp * rhs.imag_mp) / denominator;
    return ComplexGMP(newreal_mp, newimag_mp);
  }

  /// @brief type casting to complex double
  /// @return complex double
  operator std::complex<double>() {
    return std::complex<double>(real_mp.get_d(), imag_mp.get_d());
  }

  /// @brief  overload the abs function
  /// @return ComplexGMP
  mpf_class abs() const {
    mpf_class re_square = real_mp * real_mp;
    mpf_class im_square = imag_mp * imag_mp;
    mpf_class rho_square = re_square + im_square;

    mpf_class result;
    mpf_sqrt(result.get_mpf_t(), rho_square.get_mpf_t());

    return result;
  }
};

/// @brief Computes the recurrence coefficients from Thiele's continued fraction
///        This routine uses tabulation in order to efficienly compute the
///        matrix elements g_func(:,:)
/// @param x array of the frequencies
/// @param y array of the Wmn matrix elements
/// @param g_func[inout] recurrence matrix used to compute the parameters a_n
/// @param n number of parameters minus one (C-based indexing)
void thiele_pade_gcoeff_mp(const std::vector<ComplexGMP> &x,
                           const std::vector<ComplexGMP> &y,
                           std::vector<std::vector<ComplexGMP>> &g_func,
                           int n) {
  g_func[n][0] = y[n];
  if (n == 0)
    return;

  for (int idx = 1; idx <= n; idx++) {
    g_func[n][idx] = (g_func[idx - 1][idx - 1] - g_func[n][idx - 1]) /
                     ((x[n] - x[idx - 1]) * g_func[n][idx - 1]);
  }
}

/// @brief Gets the value of the Wmn matrices using the previously computed Pade
/// approximant
///        Here we only implement the Pade approximant evaluation
/// @param n_par number of parameters
/// @param x_ref array of the reference points
/// @param x the point to evaluate
/// @param a_par array of the input parameters
/// @return the value of the interpolant at x
template <typename any_complex>
any_complex evaluate_thiele_pade(int n_par,
                                 const std::vector<any_complex> &x_ref,
                                 const any_complex &x,
                                 const std::vector<any_complex> &a_par) {

  // Define constants
  any_complex c_one(std::complex<double>(1.0, 0.0));
  any_complex c_zero(std::complex<double>(0.0, 0.0));

  // Wallis' method coefficients
  std::vector<any_complex> acoef(n_par+1), bcoef(n_par+1);

  // Compute continues fraction
  acoef[0] = c_zero;
  acoef[1] = a_par[0];
  std::fill_n(bcoef.begin(), 2, c_one);

  for (int i_par = 0; i_par < n_par - 1 ; i_par++) {
    acoef[i_par + 2] =
        acoef[i_par + 1] + (x - x_ref[i_par]) * a_par[i_par + 1] * acoef[i_par];
    bcoef[i_par + 2] =
        bcoef[i_par + 1] + (x - x_ref[i_par]) * a_par[i_par + 1] * bcoef[i_par];
  }

  return acoef[n_par] / bcoef[n_par];
}

/// @brief Gets the Pade approximant of a meromorphic function F
/// @param n_par order of the interpolant
/// @param x_ref array of the reference points
/// @param y_ref array of the reference function values
/// @param a_par array of the interpolant parameters
/// @param do_greedy whether to use the default greedy algorithm or the naive
/// one
void thiele_pade_mp(int n_par, std::vector<ComplexGMP> &x_ref_mp,
                    const std::vector<ComplexGMP> &y_ref_mp,
                    std::vector<ComplexGMP> &a_par_mp,
                    std::optional<bool> do_greedy = std::nullopt) {
  // set floating point precision of GMP lib
  mpf_set_default_prec(PREC);

  // auxiliary variables
  bool local_do_greedy = true;
  int kdx, n_rem;
  std::vector<int> n_rem_idx;
  mpf_class deltap, pval;
  ComplexGMP pval_in, x_in, y_in;
  std::vector<std::vector<ComplexGMP>> g_func;
  std::vector<ComplexGMP> x(n_par), xtmp(n_par), ytmp(n_par);

  // whether to perform the refined Thiele's interpolation (default)
  if (do_greedy.has_value()) {
    local_do_greedy = do_greedy.value();
  }

  // initialize variables
  for (int i_par = 0; i_par < n_par; i_par++) {
    g_func.push_back(std::vector<ComplexGMP>());
    a_par_mp.push_back(ComplexGMP());
    n_rem_idx.push_back(i_par);
    for (int j_par = 0; j_par < n_par; j_par++) {
      g_func[i_par].push_back(ComplexGMP()); // init to 0 + 0i
    }
  }

  if (local_do_greedy) {
    x = x_ref_mp;
    x_ref_mp.assign(x_ref_mp.size(), ComplexGMP());

    // Finding the index that maximizes |F|
    auto it = std::max_element(
        y_ref_mp.begin(), y_ref_mp.end(),
        [](const ComplexGMP &a, const ComplexGMP &b) { return a.abs() < b.abs(); });
    kdx = std::distance(y_ref_mp.begin(), it);

    // Update indexes of non-visited points
    n_rem_idx.erase(n_rem_idx.begin() + kdx);
    n_rem = n_par - 1;

    // Add the winning point and compute generating function
    xtmp[0] = x[kdx];
    ytmp[0] = y_ref_mp[kdx];
    x_ref_mp[0] = x[kdx];

    thiele_pade_gcoeff_mp(xtmp, ytmp, g_func, 0);
    a_par_mp[0] = g_func[0][0];

    // Add remaining points ensuring min |P_i(x_{1+1}) - F(x_{i+1})|
    for (int idx = 1; idx < n_par; ++idx) {
      pval = mpf_class("1e9999"); // Huge value

      for (int jdx = 0; jdx < n_rem; ++jdx) {
        // Compute next convergent P_i(x_{i+1})
        pval_in = evaluate_thiele_pade<ComplexGMP>(idx, xtmp, x[n_rem_idx[jdx]], a_par_mp);

        deltap = (pval_in - y_ref_mp[n_rem_idx[jdx]]).abs();
        if (deltap < pval) {
          pval = deltap;
          x_in = x[n_rem_idx[jdx]];
          y_in = y_ref_mp[n_rem_idx[jdx]];
          kdx = jdx;
        }
      }

      // Update indexes of non-visited points
      n_rem_idx.erase(n_rem_idx.begin() + kdx);
      n_rem -= 1;

      // Add the winning point and recompute generating function
      x_ref_mp[idx] = x_in;
      xtmp[idx] = x_in;
      ytmp[idx] = y_in;
      thiele_pade_gcoeff_mp(xtmp, ytmp, g_func, idx);

      // Unpack parameters a_i = g_i(w_i)
      a_par_mp[idx] = g_func[idx][idx];
    }
  } else {
    for (int i_par = 0; i_par < n_par; i_par++) {
      thiele_pade_gcoeff_mp(x_ref_mp, y_ref_mp, g_func, i_par);
      a_par_mp[i_par] = g_func[i_par][i_par];
    }
  }
}

/// @brief function to compute Thiele-Pade approximations of a meromorphic
/// function
/// @param n_par order of the interpolant
/// @param x_ref array of the reference points
/// @param y_ref  array of the reference function values
/// @param x_query array of points where the function needs to be evaluated
/// @param y_query[out] array of the interpolated values at x_query
/// @param num_query number of query points
void thiele_pade_mp_api(int n_par, const std::complex<double> *x_ref,
                        const std::complex<double> *y_ref,
                        const std::complex<double> *x_query,
                        std::complex<double> *y_query, int num_query) {
  // set floating point precision of GMP lib
  mpf_set_default_prec(PREC);

  // initialize arbitrary precision arrays
  std::vector<ComplexGMP> x_ref_mp, y_ref_mp, x_query_mp, a_par_mp;
  for (int i_par = 0; i_par < n_par; i_par++) {
    x_ref_mp.push_back(x_ref[i_par]);
    y_ref_mp.push_back(y_ref[i_par]);
    a_par_mp.push_back(ComplexGMP());
  }
  for (int i_query = 0; i_query < num_query; i_query++) {
    x_query_mp.push_back(x_query[i_query]);
  }

  // Compute the coefficients a_par_mp
  thiele_pade_mp(n_par, x_ref_mp, y_ref_mp, a_par_mp);

  // Evaluate the Thiele-Pade approximation at the query points
  for (int i_query = 0; i_query < num_query; i_query++) {
    auto y = evaluate_thiele_pade<ComplexGMP>(n_par, x_ref_mp,
                                              x_query_mp[i_query], a_par_mp);
    y_query[i_query] = (std::complex<double>)y;
  }
}
