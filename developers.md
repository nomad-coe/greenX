# For Developers

## Portability

* The minimax library has been tested on:
  - macOS Catalina, with GCC9 and openBLAS
  - Ubuntu-20.04, with GCC9 and blas/lapack
  
* There has been no testing with Intel or MKL. This should be done.

## Testing 

The test framework is currently set up such that one:
* Writes a .f90 test driver, which calls some GreenX API
* Writes a .py file to: 
    * Supply inputs to the .f90 driver
    * Execute the driver binary
    * Parse the outputs of the test
    * Make assertions
    
This has some requirements, and some technical open questions associated with it.
Two better approaches would be to:  
a). Replace a fortran unit testing framework. In most cases, this is more appropriate.
     This would remove the need for all I/O and path-setting.
b). Call the respective fortran library directly from python, removing the need
    for the .f90 driver, and removing all the complications associated with
    consistent binary/file paths.

### Issues with the Test Set Up

#### Test Binary Name

* The fortran binary name is specified in the python test.
  - If one enforces the convention that the .f90 and .py files must have the same
    prefix, the binary name can automatically be retrieved by python
    `os.path.basename(__file__).py`

#### Location of a Test Binary

* The build directory root defined by python, from which to start searching for the 
  test binary is not robust w.r.t change in directory nesting. As directory nesting 
  is changed, one is forced to modify `root = Path(__file__).parent.parent` to 
  `root = Path(__file__).parent.parent.parent`, etc.
*`find_test_binary` recursively searches for the binary, per test!
   This will get inefficient very quickly, and should be viewed as a *place holder*.
   - A prior implementation passed the full binary path to pytest, however
   one must configure `conftest.py` such that it does not interfere with pytest's
   existing command-line arguments. This is more robust than the current 
   implementation. 
      
* Alex *assumes* relative paths will be a function of the directory in which pytest
  is invoked. If one wants to to retain the flexbility of executing pytest in 
  `<BUILD_DIR>`, `<BUILD_DIR>/test` or lower, relative binary paths are not a
  good option. 
  - One should investigate placing test binaries next to the python
    test drivers, in the build folder. Perhaps relative paths can be used, 
    removing the need to define `root` and `find_test_binary`. 

* One shuold test the library by directly calling from python, however Alex 
  (thinks) memory allocation in the library prevents this.
  - Refactor the library, such that the caller handles the memory. 
