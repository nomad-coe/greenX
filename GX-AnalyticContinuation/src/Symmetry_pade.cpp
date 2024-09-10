
// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include "Symmetry_pade.hpp"

/// @brief contruct symmetry from integer symmetry label
/// @param symm symmetry label
Symmetry_pade::Symmetry_pade(int symm){
    if ((symm >= 0) && (symm <= 7)) {
        my_symmetry = symm;
    } else {
        throw std::invalid_argument("symmetry specifier out of range!");
    }
}

/// @brief setting symmetry
/// @param symm symmetry label
void Symmetry_pade::set_symmetry(int symm){
    if ((symm >= 0) && (symm <= 7)) {
        my_symmetry = symm;
    } else {
        throw std::invalid_argument("symmetry specifier out of range!");
    }
}

/// @brief enforce pade symmetry by projecting the function argument
/// @param x function argument 
/// @return projected function argument
std::complex<double> Symmetry_pade::apply_sym_x(std::complex<double> x) {
    std::complex<double> x_symm;
    switch (my_symmetry) {
        case sym_none: {
            x_symm = x;
            break;
        }
        case sym_y: {
            x_symm = std::complex<double>(abs(x.real()), x.imag());
            break;
        }
        case sym_x: {
            x_symm = std::complex<double>(x.real(), abs(x.imag()));
            break;
        }
        case sym_xy: {
            x_symm = std::complex<double>(abs(x.real()), abs(x.imag()));
            break;
        }
        case sym_even: {
            double new_im = std::copysign(1.0, x.real()) * x.imag();
            x_symm = std::complex<double>(abs(x.real()), new_im);
            break;
        }
        case sym_odd: {
            double new_im = std::copysign(1.0, x.real()) * x.imag();
            x_symm = std::complex<double>(abs(x.real()), new_im);
            break;
        }
        case sym_conjugate: {
            double new_im = std::copysign(1.0, x.real()) * x.imag();
            x_symm = std::complex<double>(abs(x.real()), new_im);
            break;
        }
        case sym_anti_conjugate: {
            double new_im = std::copysign(1.0, x.real()) * x.imag();
            x_symm = std::complex<double>(abs(x.real()), new_im);
            break;
        }
    }
    return x_symm;
}

/// @brief enforce pade symmetry by projecting the function value if function 
///        argument was projected as well
/// @param x_original function argument before projection
/// @param y function value
/// @return projected function value
std::complex<double> Symmetry_pade::apply_sym_y(std::complex<double> x_original, std::complex<double> y) {
    std::complex<double> x_projected = apply_sym_x(x_original);
    std::complex<double> y_symm;
    if (x_projected == x_original) {
        y_symm = y;
    } else {
        switch (my_symmetry) {
            case sym_odd: {
                y_symm = -y;
                break;
            }
            case sym_conjugate: {
                y_symm = std::conj(y);
                break;
            }
            case sym_anti_conjugate: {
                y_symm = -std::conj(y);
                break;
            }
            default: {
                y_symm = y;
            }
        }
    }
    return y_symm;
}