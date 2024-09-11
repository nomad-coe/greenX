// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include "Symmetry_pade.hpp"
#include "ComplexGMP.hpp"
#include <gmpxx.h>

#include <complex>

/// @brief stores all pade parameters
struct pade_model {
    /// number of pade parameters
    int n_par;
    /// floating point precision in bit
    int precision;
    /// symmetry of pade interpolant
    Symmetry_pade symmetry;
    /// pade parameters
    ComplexGMP *a_par;
    /// reference points (resorted if greedy was used)
    ComplexGMP *xref;
};

// function prototype for fortran interface
extern "C" {

std::complex<double> evaluate_thiele_pade_mp(const std::complex<double> x,
                                             pade_model *params_ptr);

pade_model *thiele_pade_mp(int n_par, const std::complex<double> *x_ref,
                           const std::complex<double> *y_ref, int do_greedy, 
                           int precision, int symmetry);

void free_pade_model(pade_model *model) {
    if (model != nullptr) {
        delete[] model->a_par;
        delete[] model->xref;
        delete model;
    }
}

}