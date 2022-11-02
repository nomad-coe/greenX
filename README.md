# GreenX Library 

The Green X library is developed under Work Package 2 of the NOMAD Center of Excellence. 

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

If all requirements are found, build and install the project:

 ```bash
make -j
make install 
 ```

## Running the Tests

GreenX uses pytest as its regression testing framework, in conjunction with 
some custom python modules. First, one must ensure that the python utilities
are installed. From the GreenX root directory:

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

GreenX is documented using Doxygen, and documentation support is enabled by
default. To disable CMake looking for Doxygen, configure with:

```bash
cmake ../ -DENABLE_DOCS=OFF
```

To build the document, type in the build directory:

```bash
make docs
```

Documentation is built in `documentation` and can be viewed by opening
`html/index.html` in a browser.

When adding new files with documentation, please ensure the directory is listed 
in the `INPUT` tag of Doxyfile.
