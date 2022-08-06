# GreenX Library 

The Green X library is developed under Work Package 2 of the NOMAD Center of Excellence. 

## Installation

Green X has been designed as a collection of libraries, which can be built relatively 
independently of one another. To build the whole suite of Green X libraries from the source 
you need to have a Fortran compiler supporting Fortran 2008, and one of the supported build 
systems:

* cmake version 3.0.2 or newer, with a build-system backend, i.e. `make`.

**Building with CMake**   

To build all libraries, set up a build directory:

```bash
mkdir build 
```

Change to the directory and run cmake configuration:

```bash
cmake ../
```

If all requirements are found, build and install the project:

 ```bash
make -j
make install 
 ```

**Running the Tests** 

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


# TODOs 
## Issues with Robustness
* Fortran binary name specified in the python test
* root specification in python driver is not robust w.r.t change in 
  directory nesting. `root = Path(__file__).parent.parent.parent`
* `find_test_binary` recursively searches for the binary, per test!
   This will get inefficient very quickly, and should be viewed as a place holder.
   One prior implementation passed the full binary path to pytest, however
   one must configure conftest such that it does not interfere with pytest's
   existing command-line arguments. 
* I assume relative binary paths will be a function of which directory pytest
  is executed from, making them not an option (want to retain the flexbility
  to execute tests in <BUILD_DIR>, <BUILD_DIR>/test or lower). 
* Would like to do away with python writing to file, fortran reading parameters
from file.The obvious approach is to move to a unit test framework 
  In principle, one should not need to use pytest for this project. 

* Make sure shared libraries can also be built

* In time-frequency CMake, turn the application test into a function, so
  it's compact and reusable.
* Note GX lib entry point
* Would like to test the lib by directly calling fortran from python,
  however I (think) memory allocation of returned arrays by the lib prevents
  this(?) - look more into cython ptr usage.
* Document how to add new tests. 
* Check how to build the libs separately
* Check the installation works

* Update the github actions (if they've broken)
