// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include <gmpxx.h>

#include <complex>

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
    ComplexGMP(std::complex<double> z)
        : real_mp(real(z)), imag_mp(imag(z)) {}

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
    ComplexGMP operator*(const ComplexGMP &rhs) const;

    /// @brief division of two arbitrary precision complex numbers
    /// @param rhs denominator
    /// @return arbitrary precision comoplex number
    ComplexGMP operator/(const ComplexGMP &rhs) const; 

    /// @brief type casting to complex double
    /// @return complex double
    operator std::complex<double>() const {
        return std::complex<double>(real_mp.get_d(), imag_mp.get_d());
    }

    /// @brief  overload the abs function
    /// @return returns the absolut as GMP float
    mpf_class abs() const; 
};