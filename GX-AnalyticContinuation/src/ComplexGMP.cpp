// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include "ComplexGMP.hpp"

/// @brief multiplication of two arbitrary precision complex numbers
/// @param rhs right factor
/// @return arbitrary precision comoplex number
ComplexGMP ComplexGMP::operator*(const ComplexGMP &rhs) const {
    mpf_class newreal_mp = real_mp * rhs.real_mp - imag_mp * rhs.imag_mp;
    mpf_class newimag_mp = real_mp * rhs.imag_mp + imag_mp * rhs.real_mp;
    return ComplexGMP(newreal_mp, newimag_mp);
}

/// @brief division of two arbitrary precision complex numbers
/// @param rhs denominator
/// @return arbitrary precision comoplex number
ComplexGMP ComplexGMP::operator/(const ComplexGMP &rhs) const {
    mpf_class denominator =
        rhs.real_mp * rhs.real_mp + rhs.imag_mp * rhs.imag_mp;
    mpf_class newreal_mp =
        (real_mp * rhs.real_mp + rhs.imag_mp * imag_mp) / denominator;
    mpf_class newimag_mp =
        (rhs.real_mp * imag_mp - real_mp * rhs.imag_mp) / denominator;
    return ComplexGMP(newreal_mp, newimag_mp);
}


/// @brief  overload the abs function
/// @return returns the absolut as GMP float
mpf_class ComplexGMP::abs() const {
    mpf_class re_square = real_mp * real_mp;
    mpf_class im_square = imag_mp * imag_mp;
    mpf_class rho_square = re_square + im_square;
    mpf_class result;
    mpf_sqrt(result.get_mpf_t(), rho_square.get_mpf_t());
    return result;
}