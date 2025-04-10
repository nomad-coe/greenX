# GreenX library - Analytic Continuation 

The analytic continuation component (GX-AC) provides routines to interpolate functions using the thiele pade interpolant.

> [!Note]
> **Key Features**
> - basic thiele pade algorithm
> - greedy algorithm for thiele pade to enhance numerical stability
> - arbitrary precision arithmetic using the GMP library for even more numerical stability

## Building

If you want to compile only the Analytic Continuation (AC) component of Greenx, change to the GreenX root, then type:
```bash
mkdir build && cd build 
cmake -DMINIMAX_COMPONENT=OFF ..
make -j 
make install 

# optional: for running the regression tests
ctest
```

### Linking against GNU Multiple Precision (GMP) Library

If requested, the arithmetic operations to obtain the pade model are carried out using an user specified precision. These arbitrary precision floats are handeled by the [GMP library](https://gmplib.org/). **By default GreenX tries to find and link GMP automatically** since it is installed already on most systems that use GNU compilers. However, if GreenX is unable to find the GMP library it will note the user during the cmake configuration step. 

#### Compile and link GMP manually:
If for some reason no GMP is found in your system you can build it manually and specify its path in the GreenX build. Assumes that `$GX_ROOT` is set to the root of GreenX:
```bash
cd $GX_ROOT

# ------------------- Build GMP library
mkdir external && cd external 
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz 
cd gmp-6.3.0/ && mkdir install
GMP_INSTALL_DIR="$(realpath install)"
./configure --prefix="${GMP_INSTALL_DIR}" --enable-cxx
make -j 
make install 
# optionally run unit tests of GMP
#   make check

cd $GX_ROOT

# ------------------- Build GreenX
mkdir build && cd build 
cmake -DCMAKE_PREFIX_PATH=$GMP_INSTALL_DIR ../
make -j 
make install
```


#### Turn off GMP:
You can turn of linking to GMP by specifiying `ENABLE_GNU_GMP=OFF` in the cmake configuration step.
```bash
cmake -DENABLE_GNU_GMP=OFF ../              # Default: ENABLE_GNU_GMP=ON
```

## Calling the Analytic Continuation Component  

Information on how to use this component is always up-to-date on the [GreenX website](https://nomad-coe.github.io/greenX/gx_ac.html). 

Additionally, take a look at the example script in `GX-AnalyticContinuation/examples`. This script shows how this component can be used. 

## Build the Example Script 
In the `GX-AnalyticContinuation/examples` folder you can find a stand-alone fortran program (`pade_example.f90`) showcasing the usage of the GX-AC component. After building the GreenX library, change into the `GX-AnalyticContinuation/examples` folder and type:
```bash
make 
```
The program can be executed by running:
```bash 
./pade_example > output.dat
```
To plot the reference function and the Padé interpolation call the python script:
```bash 
python plot.py output.dat comparison.png
```
The figure is saved as `comparison.png`. Feel free to change some parameters in the `pade_example.f90` script and compile again to see how it is affecting the pade interpolation.

## Running the regression tests 
Regression tests of the GX-AC component use the testing framework [pytest](https://docs.pytest.org/en/stable/#). Simply type `ctest` in the build directory after the GreenX library has been build. Note that if you link against the GMP libaray (`-DENABLE_GNU_GMP=ON`), additional regression tests will run that test the behaviour of the AC component in combination with GMP.

For more information please refer to the main [README.md](https://github.com/nomad-coe/greenX/blob/main/README.md) of this repository.


