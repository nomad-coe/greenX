## GreenX library - LocalizedBasis 

This library provides a list of localized basis set procedures, commonly occurring in post-scf calculations like MP2, RPA and GW approaches. 

At the current stage, localized basis component of the GreenX library provide the separable resolution of the identity fitting coefficients to be used for the calculation of exact exchange energy, canonical RPA and GW methodologies. Also it can be used for the cubic-scaling RPA and GW, this feature is under development.

See the [GreenX document](Documents/GreenX.md) for the structure of the library. 

## Building

With CMake, change to the GreenX root, then type:

```bash
mkdir build && cd build
cmake ../
make -j
make
```


## Running the Tests

Application tests are run with the pytest framework. Having installed `pygreenx`
(see the top-level [README](../README.md)) change to `<GX_ROOT>/<BUILD_DIR>` and type `ctest`.

## Calling the Localize basis Library

Information on how to use the library is always up-to-date on the [GreenX website](nomad-coe.github.io/greenX/).

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
switched to using a fixture defined in `GX-LocalizedBasis/test/conftest.py`,
for the duration of the test execution. The `tmp_path` fixture will provide a
temporary directory unique to the (py)test invocation, created in the base
temporary directory.

A new test can be added to this `CMakeLists.txt` in the following way:

```cmake
set(LIBS_FOR_TESTING)
list(APPEND LIBS_FOR_TESTING LibGXLBasis)

set(TEST_NAME "test_gx_localize_basis")
add_app_test(TEST_NAME TEST_TARGET_DIR LIBS_FOR_TESTING)
```

Calling `add_app_test` should ensure that the directory level structure in
the build folder is consistent with example python driver. For multiple tests,
one can place the call to `add_app_test` in a loop, and supply different test names.

All the tests can be run by switching to the `build/` directory and typing `ctest`.
