#include <complex>
#include <gmpxx.h>

// Set GMP precision
const int PREC = 128;

/// @brief  datatype of a complex number using arbitrary precision numbers for
///         real and imaginary part
class ComplexGMP {
public:
  /// real part
  mpf_class real_mp;
  /// imag part
  mpf_class imag_mp;

  /// @brief initialize a complex number 0 + 0i
  ComplexGMP() : real_mp(0, PREC), imag_mp(0, PREC) {}

  /// @brief initialize a complex number from real and imag part using arbitrary
  ///        precision numbers
  /// @param real_mp real part arbitrary precision
  /// @param imag_mp imaginary part arbitrary precision
  ComplexGMP(const mpf_class &real_mp, const mpf_class &imag_mp)
      : real_mp(real_mp, PREC), imag_mp(imag_mp, PREC) {}

  /// @brief initialize a complex number from a complex double
  /// @param z complex number
  ComplexGMP(std::complex<double> z)
      : real_mp(real(z), PREC), imag_mp(imag(z), PREC) {}

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

struct pade_model {
  int n_par;
  ComplexGMP *a_par;
  ComplexGMP *xref;
};

// function prototype for fortran interface
extern "C" {
std::complex<double> evaluate_thiele_pade_mp(const std::complex<double> x,
                                             pade_model *params_ptr);
pade_model *thiele_pade_mp(int n_par, const std::complex<double> *x_ref,
                           const std::complex<double> *y_ref, int do_greedy);
void free_pade_model(pade_model *model) {
  if (model != nullptr) {
    delete[] model->a_par;
    delete[] model->xref;
    delete model;
  }
}
}