#include <complex>
#include <gmpxx.h>

class ComplexGMP;

struct pade_model {
    int n_par;
    ComplexGMP *a_par;
    ComplexGMP *xref;
};

// function prototype for fortran interface
extern "C" {
std::complex<double> evaluate_thiele_pade_mp(const std::complex<double> x,
                                             pade_model *params_ptr);
pade_model *thiele_pade_mp(int n_par, 
                           const std::complex<double> *x_ref,
                           const std::complex<double> *y_ref, 
                           int do_greedy);
}

