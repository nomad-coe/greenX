<h1 align="center">
  <img src="docs/logo/greenx_logo.png" alt="GreenX" width="300">
</h1>

# GreenX Library 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.05570/status.svg)](https://doi.org/10.21105/joss.05570)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07859/status.svg)](https://doi.org/10.21105/joss.07859)

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

Specific shared libraries can be disabled or enabled, e.g to enable the Projector-Augmented Wave (PAW) component of GreenX, run CMake configure with: 

```bash
cmake ../ -DPAW_COMPONENT=ON
```

Available options to disable one or more components

| Component                              | CMake configure                    | Default |
|----------------------------------------|------------------------------------|---------|
| Analytical Continuation component      | `-DAC_COMPONENT`      | `ON`  |
| Minimax Time-Frequency grids component | `-DMINIMAX_COMPONENT` | `ON`  |
| Localized Basis component              | `-DLBASIS_COMPONENT`  | `OFF` |
| Projector-Augmented Wave component     | `-DPAW_COMPONENT`     | `OFF` |

GreenX uses GNU Multiple Precision Arithmetic Library by default in the Analytical Continuation component, you can disable it without any harm by running CMake configure with:

```bash
cmake ../ -DENABLE_GNU_GMP=OFF
```

To build GreenX with submodules (they are turned off by default), one can configure with:

```bash
cmake ../ -DCOMPILE_SUBMODULES=ON
```

You can obtained them by executing:

```bash
git submodule update --init --recursive
```

If all requirements are found, build and install the project:

 ```bash
make -j
make install 
 ```

Minimal example to build and install only the Time-Frequency component of GreenX:

 ```bash
cmake .. -DAC_COMPONENT=OFF 
make -j
make install
 ```


## Regression Tests

### Running Regression Tests
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

#### Adding Regression Tests 


 1. Add a New Test Source File: Create a test file in the `test/` directory of the component. If testing a Fortran function, create a file like `test_new_feature.f90` with the necessary test logic.

 2. Update `CMakeLists.txt`: Modify `CMakeLists.txt` in the component directory to include the new test:
    - Add the source file to the test executable
        ```cmake
        # Define the new test target
        add_executable(test_gx_new_feature)

        # Set binary name
        set_target_properties(test_gx_new_feature
            PROPERTIES
            RUNTIME_OUTPUT_NAME test_gx_new_feature.exe)

        # Add source file for the new test
        target_sources(test_gx_new_feature
            PRIVATE
            test/test_new_feature.f90
        )

        # Link the test executable to the necessary libraries
        target_link_libraries(test_gx_new_feature
            PUBLIC
            LibGXNewFeature
        )

        # Specify the runtime output directory
        set_target_properties(test_gx_new_feature
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})
        ```

    - Copy Python test files. If the test involves Python, ensure the test scripts are copied to the build directory:
        ```cmake
        add_custom_command(
            TARGET LibGXNewFeature POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_new_feature.py
                    ${PROJECT_BINARY_DIR}/test/new_feature/)
        ```

    - Add the new test to CTest
        ```cmake
        add_test(
            NAME test_gx_new_feature
            COMMAND pytest -s test_new_feature.py --build-dir ${CMAKE_BINARY_DIR}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/new_feature
        )
        ```

### Running Unit Tests

The GX-q=0 component of GreenX uses the unit-testing framework [Zofu](https://github.com/acroucher/zofu). This library is build together with Greenx when `ENABLE_GREENX_UNIT_TESTS=ON`: 

```bash 
cmake -DENABLE_GREENX_UNIT_TESTS=ON ../
```
Unit tests are run with the application tests, using ctest. Simply type `ctest`
in the build directory.

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


## Contribute

Contributions to GreenX are highly appreciated! If you consider contributing [check out CONTRIBUTING.md](https://github.com/nomad-coe/greenX/blob/main/CONTRIBUTING.md).
