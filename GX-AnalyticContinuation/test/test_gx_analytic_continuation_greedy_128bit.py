"""
**************************************************************************************************
 Copyright (C) 2020-2024 GreenX library
 This file is distributed under the terms of the APACHE2 License.

**************************************************************************************************

Asertions for test_gx_analytic_continuation.f90

"""

import numpy as np
from pathlib import Path
import pytest
import os

from pygreenx.run import BinaryRunner, BuildType

@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'test_gx_analytic_continuation.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary


def get_AC_residual_sum(fortran_binary):
    """ run the AC component"""
    # Run test
    runner = BinaryRunner(fortran_binary, BuildType.serial, args=["greedy", "128", "none", "normal"])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"

    # get result
    residual_sum = np.loadtxt("residual_sum.txt", comments="#")

    return residual_sum


def test_get_rirs_overlap_error(fortran_binary): 
    """ test the AC component """
    ref =  0.61478063E-06                         # reference value
    res_sum = get_AC_residual_sum(fortran_binary) # comupte actual value
    print(f"comparing ref:{ref}, out:{res_sum}")
    assert np.allclose(res_sum, ref, atol=1.e-8)

