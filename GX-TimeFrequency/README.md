# GreenX Library - TimeFrequency

## Building

With CMake, change to the GreenX root, then type:

```bash
mkdir build && cd build
cmake ../
make -j 
make install 
```

## Running the Tests

Application tests are run with the pytest framework. Having installed `pygreenx`
(see the top-level README) change to `<GX_ROOT>/<BUILD_DIR>` and type `ctest`. 
Additionally, one can change to the test folder and explicitly run the pytest 
command from there:

```bash
cd <GX_ROOT>/<BUILD_DIR>/test/time-frequency
pytest -s 
# or
pytest -s test_name.py
```

## Calling the MiniMax Library

The minimax grid generation is called like so:

```fortran
use gx_minimax, only: gx_minimax_grid

! Declarations
integer :: n_mesh_points
real(dp) :: e_transition_min, e_transition_max
real(dp), allocatable :: tau_mesh(:), tau_weights(:)
real(dp), allocatable :: freq_mesh(:), freq_weights(:)
real(dp), allocatable :: cos_tau_to_freq_weights(:, :)
real(dp), allocatable :: cos_freq_to_tau_weights(:, :)
real(dp), allocatable :: sinft_tau_to_freq_weights(:, :)
real(dp) :: max_errors(3)
real(dp) :: cosft_duality_error
integer :: ierr

call gx_minimax_grid(n_mesh_points, e_transition_min, e_transition_max, &
                     tau_mesh, tau_weights, &
                     freq_mesh, freq_weights, &
                     cos_tau_to_freq_weights, cos_freq_to_tau_weights, &
                     sinft_tau_to_freq_weights, &
                     max_errors, cosft_duality_error, ierr)
```

For a description of the variables, please consult `src/minimax_grids.F90`.
For an example of how to call `gx_minimax_grid`, please consult `test/test_gx_minimax_grid.f90`.

Additionally, one can also call a utility routine to query whether the 
number of imaginary-time points has a corresponding grid tabulation:

```fortran
use api_utilites, only: gx_check_ntau, gx_get_error_message

call gx_check_ntau(ntau, msg, ierr)
```

## For Developers

### Adding a New Test

Write a `test_name.py` and `test_name.f90` in `test/`. 
`test_name.f90` should accept some input parameter/s, make a library call and
write the result to file. `test_name.py` should provide the input parameters
(the simplest way is to write to file), execute the binary built from 
`test_name.f90`, parse the result and make the assertions. This implies the
reference results should also be stored in `test_name.py`.

`test_name.py` *must* include the fixture:

```python
@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'gx_tabulate_grids.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary
```

where the binary `name` should be consistent with whatever the developer has
called the binary in the CMakeLists.txt

One notes that the working directory for all tests is `tmp_path`, which is 
switched to using a fixture defined in `GX-TimeFrequency/test/conftest.py`, 
for the duration of the test execution. The `tmp_path` fixture will provide a 
temporary directory unique to the (py)test invocation, created in the base 
temporary directory.

A new test can be added to this `CMakeLists.txt` in the following way:

```cmake
set(LIBS_FOR_TESTING)
list(APPEND LIBS_FOR_TESTING LibGXMiniMax)

set(TEST_NAME "test_gx_minimax_grid")
add_app_test(TEST_NAME TEST_TARGET_DIR LIBS_FOR_TESTING)
```

Calling `add_app_test` should ensure that the directory level structure in
the build folder is consistent with example python driver. For multiple tests, 
one can place the call to `add_app_test` in a loop, and supply different test names.

All the tests can be run by switching to the `build/` directory and typing `ctest`,
however they can also be run in the `test/` directory of the source by typing
`pytest -s` as long as the 'ENV VAR' to the build directory has been defined:

```bash
export GX_BUILD_DIR=<path/2/build/>
```
