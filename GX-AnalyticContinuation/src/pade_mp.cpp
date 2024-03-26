// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include "pade_mp.h"

#include <complex>
#include <iostream>
#include <optional>
#include <vector>

/// @brief Computes the recurrence coefficients from Thiele's continued fraction
///        This routine uses tabulation in order to efficienly compute the
///        matrix elements g_func(:,:)
/// @param x array of the frequencies
/// @param y array of the Wmn matrix elements
/// @param g_func[inout] recurrence matrix used to compute the parameters a_n
/// @param n number of parameters minus one (C-based indexing)
void thiele_pade_gcoeff_mp(const ComplexGMP *x,
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
///        approximant
///        Here we only implement the Pade approximant evaluation
/// @param n_par number of parameters
/// @param x_ref array of the reference points
/// @param x the point to evaluate
/// @param a_par array of the input parameters
/// @return the value of the interpolant at x
template <typename any_complex>
void evaluate_thiele_pade_tab_mp(int n_par, const any_complex *x_ref,
                                 const any_complex &x, const any_complex *a_par,
                                 std::vector<any_complex> &acoef,
                                 std::vector<any_complex> &bcoef) {
    // Define constants and variables
    const ComplexGMP c_one(std::complex<double>(1.0, 0.0));
    any_complex delta;

    // Wallis' method iteration
    delta = a_par[n_par + 1] * (x - x_ref[n_par]);
    if (n_par == 0) {
        acoef[1] = acoef[0];
        bcoef[1] = c_one + delta;
        return;
    }
    acoef[n_par + 1] = acoef[n_par] + delta * acoef[n_par - 1];
    bcoef[n_par + 1] = bcoef[n_par] + delta * bcoef[n_par - 1];
}

/// @brief Evaluates a Pade approximant constructed with Thiele's reciprocal
/// differences
/// @param x the point to evaluate
/// @param params_ptr pointer to the struct holding model parameters
/// @return the value of the interpolant at x
std::complex<double> evaluate_thiele_pade_mp(const std::complex<double> x,
                                             pade_model *params) {
    mpf_set_default_prec(params->precision);

    // conversion from double to gmp
    ComplexGMP x_mp(x);

    // Define constants
    const ComplexGMP c_one(std::complex<double>(1.0, 0.0));
    const mpf_class tol("1", params->precision);

    // Wallis' method coefficients
    std::vector<ComplexGMP> acoef(params->n_par + 1), bcoef(params->n_par + 1);
    acoef[0] = params->a_par[0];
    bcoef[0] = c_one;

    // Compute continued fraction
    for (int i_par = 0; i_par < params->n_par - 1; ++i_par) {
        evaluate_thiele_pade_tab_mp<ComplexGMP>(i_par, params->xref, x_mp,
                                                params->a_par, acoef, bcoef);
        if (bcoef[i_par + 1].abs() > tol) {
            acoef[i_par + 1] = acoef[i_par + 1] / bcoef[i_par + 1];
            acoef[i_par] = acoef[i_par] / bcoef[i_par + 1];
            bcoef[i_par] = bcoef[i_par] / bcoef[i_par + 1];
            bcoef[i_par + 1] = c_one;
        }
    }

    return std::complex<double>(acoef[params->n_par - 1] /
                                bcoef[params->n_par - 1]);
}

/// @brief Gets the Pade approximant of a meromorphic function F
/// @param n_par order of the interpolant
/// @param x_ref array of the reference points
/// @param y_ref array of the reference function values
/// @param do_greedy whether to use the default greedy algorithm or the naive
/// @param precision floating point arythmetic precision in bits (!! not bytes!!)
/// @return a pointer to the struct holding all model parameters
pade_model *thiele_pade_mp(int n_par, const std::complex<double> *x_ref,
                           const std::complex<double> *y_ref, int do_greedy, int precision) {
    // set floating point precision of GMP lib
    mpf_set_default_prec(precision);

    // Define constants
    const ComplexGMP c_one(std::complex<double>(1.0, 0.0));
    const ComplexGMP c_zero(std::complex<double>(0.0, 0.0));
    const mpf_class tol("1e-6", precision);

    // auxiliary variables
    int kdx, n_rem;
    std::vector<int> n_rem_idx;
    mpf_class deltap, pval;
    ComplexGMP pval_in, x_in, y_in, acoef_in, bcoef_in;
    std::vector<std::vector<ComplexGMP>> g_func;
    ComplexGMP *a_par_mp = new ComplexGMP[n_par];
    ComplexGMP *x_ref_mp = new ComplexGMP[n_par];
    ComplexGMP *xtmp = new ComplexGMP[n_par];
    std::vector<ComplexGMP> acoef(n_par + 1), bcoef(n_par + 1);
    std::vector<ComplexGMP> y_ref_mp;
    std::vector<ComplexGMP> x(n_par), ytmp(n_par);

    // initialize variables
    for (int i_par = 0; i_par < n_par; i_par++) {
        x_ref_mp[i_par] = x_ref[i_par];
        x[i_par] = x_ref[i_par];
        y_ref_mp.push_back(y_ref[i_par]);
        g_func.push_back(std::vector<ComplexGMP>());
        a_par_mp[i_par] = ComplexGMP();
        n_rem_idx.push_back(i_par);
        for (int j_par = 0; j_par < n_par; j_par++) {
            g_func[i_par].push_back(ComplexGMP());
        }
    }

    if (do_greedy == 1) {
        for (int i_par = 0; i_par < n_par; i_par++) {
            x_ref_mp[i_par] = ComplexGMP();
        }

        /**********************************
                   a_0 coefficient
        ***********************************/

        // Finding the index that maximizes |F|
        auto it_zero =
            std::max_element(y_ref_mp.begin(), y_ref_mp.end(),
                             [](const ComplexGMP &a, const ComplexGMP &b) {
                                 return a.abs() < b.abs();
                             });
        int kdx_zero = std::distance(y_ref_mp.begin(), it_zero);

        // Add slected point
        xtmp[0] = x[kdx_zero];
        ytmp[0] = y_ref_mp[kdx_zero];
        x_ref_mp[0] = x[kdx_zero];

        // Compute the generating function
        thiele_pade_gcoeff_mp(xtmp, ytmp, g_func, 0);
        a_par_mp[0] = g_func[0][0];

        /**********************************
                   a_1 coefficient
        ***********************************/

        // Finding the index that maximizes abs(x-x_0) while excluding the first_kdx
        int kdx_one = kdx_zero + 1;
        ComplexGMP max_one;
        for (int i = 0; i < x.size(); ++i) {
            if (i == kdx_zero)
                continue;

            ComplexGMP current_diff = x[i] - x_ref_mp[0];
            if (i == 0 || current_diff.abs() > max_one.abs()) {
                kdx_one = i;
                max_one = current_diff;
            }
        }

        // Add selected point
        xtmp[1] = x[kdx_one];
        ytmp[1] = y_ref_mp[kdx_one];
        x_ref_mp[1] = x[kdx_one];

        // Compute the generating function
        thiele_pade_gcoeff_mp(xtmp, ytmp, g_func, 1);
        a_par_mp[1] = g_func[1][1];

        /**********************************
              Wallis' method coefficients
        ***********************************/

        acoef[0] = a_par_mp[0];
        bcoef[0] = c_one;

        // Update indexes of non-visited points
        if (kdx_zero > kdx_one)
            std::swap(kdx_zero, kdx_one);

        n_rem_idx.erase(n_rem_idx.begin() + kdx_one);
        n_rem_idx.erase(n_rem_idx.begin() + kdx_zero);
        n_rem = n_par - 2;

        // Add remaining points ensuring min |P_i(x_{1+1}) - F(x_{i+1})|
        for (int idx = 2; idx < n_par; ++idx) {
            pval = mpf_class("1e9999", precision);  // Huge value

            for (int jdx = 0; jdx < n_rem; ++jdx) {
                // Compute next convergent P_i(x_{i+1})
                evaluate_thiele_pade_tab_mp<ComplexGMP>(
                    idx - 2, xtmp, x[n_rem_idx[jdx]], a_par_mp, acoef, bcoef);
                pval_in = acoef[idx - 1] / bcoef[idx - 1];
                deltap = (pval_in - y_ref_mp[n_rem_idx[jdx]]).abs();

                if (deltap < pval) {
                    pval = deltap;
                    x_in = x[n_rem_idx[jdx]];
                    y_in = y_ref_mp[n_rem_idx[jdx]];
                    acoef_in = acoef[idx - 1];
                    bcoef_in = bcoef[idx - 1];
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

            // Rescale Wallis coefficients to avoid overflow
            acoef[idx - 1] = acoef_in;
            bcoef[idx - 1] = bcoef_in;
            if (bcoef_in.abs() > tol) {
                acoef[idx - 1] = acoef[idx - 1] / bcoef_in;
                acoef[idx - 2] = acoef[idx - 2] / bcoef_in;
                bcoef[idx - 2] = bcoef[idx - 2] / bcoef_in;
                bcoef[idx - 1] = c_one;
            }

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

    // Create pointer to the parameter struct
    pade_model *params = new pade_model;
    params->a_par = a_par_mp;
    params->precision = precision;
    params->n_par = n_par;
    params->xref = x_ref_mp;

    // Clean-up
    delete[] xtmp;

    return params;
}
