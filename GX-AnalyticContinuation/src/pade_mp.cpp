#include <gmpxx.h>

#include <complex>
#include <iostream>
#include <vector>

const int PREC = 1024;

// function prototype for fortran interface
extern "C" {
void thiele_pade_mp_api(int n_par,
                        const std::complex<double> *x_ref,
                        const std::complex<double> *y_ref,
                        const std::complex<double> *x_query,
                        std::complex<double> *y_query,
                        int num_query);
}

/// @brief  datatype of a complex number using arbitrary precision numbers for real and imaginary part
class ComplexGMP {
   public:
    /// real part
    mpf_class real_mp;
    /// imag part
    mpf_class imag_mp;

    /// @brief initialize a complex number 0 + 0i
    ComplexGMP() : real_mp(0), imag_mp(0) {}

    /// @brief initialize a complex number from real and imag part using arbitrary precision numbers
    /// @param real_mp real part arbitrary precision
    /// @param imag_mp imaginary part arbitrary precision
    ComplexGMP(const mpf_class &real_mp, const mpf_class &imag_mp) : real_mp(real_mp), imag_mp(imag_mp) {}

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
        mpf_class denominator = rhs.real_mp * rhs.real_mp + rhs.imag_mp * rhs.imag_mp;
        mpf_class newreal_mp = (real_mp * rhs.real_mp + rhs.imag_mp * imag_mp) / denominator;
        mpf_class newimag_mp = (rhs.real_mp * imag_mp - real_mp * rhs.imag_mp) / denominator;
        return ComplexGMP(newreal_mp, newimag_mp);
    }

    /// @brief type casting to complex double
    /// @return complex double
    operator std::complex<double>() {
        return std::complex<double>(real_mp.get_d(), imag_mp.get_d());
    }
};

/// @brief Computes the recurrence coefficients from Thiele's continued fraction.
///        This routine uses tabulation in order to efficienly compute the matrix elements g_func(:,:)
/// @param x array of the frequencies
/// @param y array of the Wmn matrix elements
/// @param g_func[inout] recurrence matrix used to compute the parameters a_n
/// @param n number of parameters minus one (C-based indexing)
void thiele_pade_gcoeff_mp(const std::vector<ComplexGMP> &x,
                           const std::vector<ComplexGMP> &y,
                           std::vector<std::vector<ComplexGMP>> &g_func,
                           int n) {
    g_func[n][0] = y[n];
    if (n == 0) return;

    for (int idx = 1; idx <= n; idx++) {
        g_func[n][idx] =
            (g_func[idx - 1][idx - 1] - g_func[n][idx - 1]) /
            ((x[n] - x[idx - 1]) * (x[n] + x[idx - 1]) * g_func[n][idx - 1]);
    }
}

/// @brief Gets the value of the Wmn matrices using the previously computed Pade approximant.
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
    any_complex c_one(std::complex<double>(1.0, 0.0));
    any_complex gtmp(c_one);

    for (int i_par = n_par - 1; i_par > 0; i_par--) {
        gtmp = c_one + a_par[i_par] * (x - x_ref[i_par - 1]) * (x + x_ref[i_par - 1]) / gtmp;
    }

    return a_par[0] / gtmp;
}

/// @brief Gets the Pade approximant of a meromorphic function F
/// @param n_par order of the interpolant
/// @param x_ref array of the reference points
/// @param y_ref array of the reference function values
/// @param a_par array of the interpolant parameters
void thiele_pade_mp(int n_par,
                    std::vector<ComplexGMP> &x_ref_mp,
                    const std::vector<ComplexGMP> &y_ref_mp,
                    std::vector<ComplexGMP> &a_par_mp) {
    // set floating point precision of GMP lib
    mpf_set_default_prec(PREC);

    // initialize array
    std::vector<std::vector<ComplexGMP>> g_func;
    for (int i_par = 0; i_par < n_par; i_par++) {
        g_func.push_back(std::vector<ComplexGMP>());
        for (int j_par = 0; j_par < n_par; j_par++) {
            g_func[i_par].push_back(ComplexGMP());  // init to 0 + 0i
        }
    }

    // interpolate
    thiele_pade_gcoeff_mp(x_ref_mp, y_ref_mp, g_func, 0);
    a_par_mp[0] = g_func[0][0];
    for (int i_par = 1; i_par < n_par; i_par++) {
        thiele_pade_gcoeff_mp(x_ref_mp, y_ref_mp, g_func, i_par);
        a_par_mp[i_par] = g_func[i_par][i_par];
    }
}

/// @brief function to compute Thiele-Pade approximations of a meromorphic function
/// @param n_par order of the interpolant
/// @param x_ref array of the reference points
/// @param y_ref  array of the reference function values
/// @param x_query array of points where the function needs to be evaluated
/// @param y_query[out] array of the interpolated values at x_query
/// @param num_query number of query points
void thiele_pade_mp_api(int n_par,
                        const std::complex<double> *x_ref,
                        const std::complex<double> *y_ref,
                        const std::complex<double> *x_query,
                        std::complex<double> *y_query,
                        int num_query) {
    // set floating point precision of GMP lib
    mpf_set_default_prec(PREC);

    // initialize arbitrary precision arrays
    std::vector<ComplexGMP> x_ref_mp, y_ref_mp, x_query_mp, a_par_mp;
    for (int i_par = 0; i_par < n_par; i_par++) {
        x_ref_mp.push_back(x_ref[i_par]);
        y_ref_mp.push_back(y_ref[i_par]);
        x_query_mp.push_back(x_query[i_par]);
        a_par_mp.push_back(ComplexGMP());
    }

    // Compute the coefficients a_par_mp
    thiele_pade_mp(n_par, x_ref_mp, y_ref_mp, a_par_mp);

    // Evaluate the Thiele-Pade approximation at the query points
    for (int i_query = 0; i_query < num_query; i_query++) {
        auto y = evaluate_thiele_pade<ComplexGMP>(n_par, x_ref_mp, x_query_mp[i_query], a_par_mp);
        y_query[i_query] = (std::complex<double>)y;
    }
}
