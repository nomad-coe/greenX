// ***************************************************************************************************
//  Copyright (C) 2020-2024 GreenX library
//  This file is distributed under the terms of the APACHE2 License.
//
// ***************************************************************************************************

#include <complex>

/// @brief utilities to enforce a special symmetry onto the pade model  
class Symmetry_pade {

    private:

        /// symmetry of the considered pade function
        int my_symmetry;

        /// symmetry type: no symmetry
        static const int sym_none = 0;
        /// symmetry type: mirror at img axis
        static const int sym_mirror_real = 1;
        /// symmetry type: mirror at real axis
        static const int sym_mirror_imag = 2;
        /// symmetry type: mirror at img and real axis
        static const int sym_mirror_both = 3;
        /// symmetry type: even function
        static const int sym_even = 4;
        /// symmetry type: odd function
        static const int sym_odd = 5;
        /// symmetry type: f(z) = \overline{f(-z)}
        static const int sym_conjugate = 6;
        /// symmetry type: f(z) = -\overline{f(-z)}
        static const int sym_anti_conjugate = 7;

    public:

        Symmetry_pade(){}
        Symmetry_pade(int symm);

        void set_symmetry(int symm);
        std::complex<double> apply_sym_x(std::complex<double> x); 
        std::complex<double> apply_sym_y(std::complex<double> x_original, std::complex<double> y); 
};