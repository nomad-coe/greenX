<h1 align="center">
  <img src="docs/logo/greenx_logo.png" alt="GreenX" width="300">
</h1>

# GreenX Library 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.05570/status.svg)](https://doi.org/10.21105/joss.05570)

The Green X library is developed under Work Package 2 of the NOMAD Center of Excellence. 
It is available under the APACHE2 [license](LICENSE.txt).

## Libraries

* [GX-AnalyticContinuation](GX-AnalyticContinuation): Performs an analytical continuation of the self-energy from the imaginary frequency to the real frequency
* [GX-LAPW](https://github.com/nomad-coe/greenX/tree/main/GX-LAPW): A cubic scaling GW algorithm in LAPW+lo basis.
* [GX-LocalizedBasis](https://github.com/nomad-coe/greenX/tree/main/GX-LocalizedBasis): The implementation of the separable resolution of the identity.
* [GX-PAW](https://github.com/nomad-coe/greenX/tree/main/GX-PAW): Supports the projector-augmented wave method 
* [GX-PlaneWaves](https://github.com/nomad-coe/greenX/tree/main/GX-PlaneWaves): Low-scaling plane-wave based GW implementation
* [GX Time-frequency](GX-TimeFrequency/README.md): Optimised quadrature grids and weights for 
  RPA and GW imaginary time-frequency transforms.
* [GX-q=0](https://github.com/nomad-coe/greenX/tree/main/GX-q%3D0): A code-agnostic framework for the treatment of inverse
dielectric function and/or screened Coulomb potential at q=0.


## Installation

Green X has been designed as a collection of libraries, which can be built relatively 
independently of one another. To build the whole suite of Green X libraries from the source 
you need to have a Fortran compiler supporting Fortran 2008, and one of the supported build 
systems:

* cmake version 3.15.0 or newer, with a build-system backend, i.e. `make`.

**Building with CMake**   

To build all libraries, set up a build directory, change to it and run cmake 
configuration:

```bash
mkdir build 
cd build
cmake ../
```

To explicitly specify the compiler, run CMake configure with:

```bash
FC=ifort cmake ../
```
Shared libraries are built by default. To build static versions, one can 
configure with:

```bash
cmake ../ -DBUILD_SHARED_LIBS=OFF
```

Specific shared libraries can be disable, e.g to disable the Projector-Augmented Wave (PAW) component of GreenX, run CMake configure with: 

```bash
cmake ../ -DPAW_COMPONENT=OFF
```

Available options to disable one or more components

| Component                              | CMake configure                    |
|----------------------------------------|------------------------------------|
| Analytical Continuation component      | `cmake .. -DAC_COMPONENT=OFF`      |
| Minimax Time-Frequency grids component | `cmake .. -DMINIMAX_COMPONENT=OFF` |
| Localized Basis component              | `cmake .. -DLBASIS_COMPONENT=OFF`  |
| Projector-Augmented Wave component     | `cmake .. -DPAW_COMPONENT=OFF`     |

GreenX uses GNU Multiple Precision Arithmetic Library by default in the Analytical Continuation component, you can disable it without any harm by runing CMake configure with:

```bash
cmake ../ -DENABLE_GNU_GMP=OFF
```

GreenX uses submodules and they are built by default. You can obtained them by executing:

```bash
git submodule update --init --recursive
```

To build GreenX without the submodules, one can configure with:

```bash
cmake ../ -DCOMPILE_SUBMODULES=OFF
```

If all requirements are found, build and install the project:

 ```bash
make -j
make install 
 ```

Minimal example to build and install only the Time-Frequency component of GreenX:

 ```bash
cmake .. -DAC_COMPONENT=OFF -DLBASIS_COMPONENT=OFF -DPAW_COMPONENT=OFF -DCOMPILE_SUBMODULES=OFF
make -j
make install
 ```


## Running the Tests

GreenX uses pytest as its regression testing framework, in conjunction with 
the custom python module `pygreenx`. First, one must ensure that `pygreenx`
is installed. From the GreenX root directory:

```bash
cd python
pip install -e .
```

One notes that the user may wish to change the scope of `pip install` to the
Python user install directory of their platform. Typically `~/.local/`. This
can be achieved with:

```bash
pip install --user -e .
```

unless working in a virtual environment, in which case it is not required.  

The test suite can now be run with CMake's ctest command:

 ```bash
cd build
ctest
 ```

## Building Documentation

GreenX is documented using Doxygen, and documentation support is disabled by
default. To enable CMake looking for Doxygen, configure with:

```bash
cmake ../ -DENABLE_GREENX_DOCS=ON
```

To build the document, type in the build directory:

```bash
make docs
```

Documentation is built in `documentation` and can be viewed by opening
`html/index.html` in a browser.

When adding new files with documentation, please ensure the directory is listed 
in the `INPUT` tag of Doxyfile.

For more information and benchmark examples see also the [GreenX website](https://nomad-coe.github.io/greenX/).

## Unit Testing

### Running Unit Tests

Unit tests require the unit-testing framework [Zofu](https://github.com/acroucher/zofu). This library is build together with Greenx when `ENABLE_GREENX_UNIT_TESTS=ON`: 

```bash 
cmake -DENABLE_GREENX_UNIT_TESTS=ON ../
```
Unit tests are run with the application tests, using ctest. Simply type `ctest`
in the build directory.


### Installing the Unit-Testing Framework manually

It is also possible to compile Zofu manually. To build Zofu, from GreenX's root (noting that one must define `$GX_ROOT`):

```bash 

# building zofu
mkdir external && cd external
git clone https://github.com/acroucher/zofu.git
cd zofu
mkdir build && cd build
cmake \
   -DCMAKE_BUILD_TYPE=release \
   -DCMAKE_INSTALL_PREFIX=${GX_ROOT}/external/zofu/install \
   -DZOFU_FORTRAN_MODULE_INSTALL_DIR:PATH=include \
   ..
make -j 4
make install

# building GreenX 
cd $GX_ROOT 
mkdir build && cd build 
cmake -DENABLE_GREENX_UNIT_TESTS=ON -DZOFU_PATH=${GX_ROOT}/external/zofu/install ../
make -j 2
make install
```
Again, typing `ctest` in the GreenX build directory starts the unit tests together with the application tests.

### Adding Unit Tests

For an example of writing a unit test, see `GX-AnalyticContinuation/src/test_pade_approximant.f90`.
A unit test is a module and follows the naming convention `test_MODULENAME.f90`.
The unit test itself is a module containing subroutines which set up some data,
call the routine under test, and make some assertions on the resulting data.
Zofu provides the object with which to make the assertions and carry the result.
One should write a separate test module for each fortran module they wish to test.

Unit tests are added to the build system straightforwardly:

1. Create a directory in the build folder that will contain the test binary.
A good convention is `unit_tests/sublibrary_name`:

```cmake
set(UNIT_TEST_DIR "${PROJECT_BINARY_DIR}/unit_tests/analytic-continuation")
file(MAKE_DIRECTORY ${UNIT_TEST_DIR})
message("-- Analytic continuation unit tests written to: ${UNIT_TEST_DIR}")
```

Noting one is free to choose any name for `UNIT_TEST_DIR`.

2. Create a list of libraries which your  unit tests depend upon. Typically
the library associated with that subfolder, for which the module is a part of.

```cmake
# Libraries on which the unit tests depend
set(LIBS_FOR_UNIT_TESTS LibGXAC)
```

Noting one is free to choose any name for `LIBS_FOR_UNIT_TESTS`.

3. Call the function `create_unit_test_executable` to define the unit test:

```cmake
create_unit_test_executable(TARGET_TEST_DIR ${UNIT_TEST_DIR}
                            TEST_NAME "test_pade_approximant"
                            REQUIRED_LIBS ${LIBS_FOR_UNIT_TESTS})
```

For multiple tests, one could call `create_unit_test_executable` in a loop over 
a list of modules.


## Contribute

Contributions to GreenX are highly appreciated! If you consider contributing [check out CONTRIBUTING.md](https://github.com/nomad-coe/greenX/blob/main/CONTRIBUTING.md).