# GreenX library - TimeFrequency

This directory is a placeholder for the development of the TimeFrequency component of GreenX.
See the [GreenX document](Documents/Green-X.md) for the structure of the library.

## Building

With CMake, change to the GreenX root, then type:

```bash
mkdir build && cd build
cmake ../
make -j 
make install 
```

## Application Tests

Application tests are run with the pytest framework. From the GreenX root directory, build the
code then type:

`pytest GX-TimeFrequency/tests/test_minimax/test_minimax.py`.
