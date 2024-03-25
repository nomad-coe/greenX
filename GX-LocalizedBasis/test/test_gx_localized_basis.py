"""
! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
Asertions for test_gx_localized_basis.f90

"""

import numpy as np
from pathlib import Path
import pytest
import os

from pygreenx.run import BinaryRunner, BuildType

@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'test_gx_localized_basis.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary


def get_rirs_overlap_error(fortran_binary):

    # Current directory

    current_dir = os.getcwd()

    # Run test
    runner = BinaryRunner(fortran_binary, BuildType.serial, args=[current_dir])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"

    # Parse results
    error_file = "error.dat"

    return np.genfromtxt(error_file)


def test_get_rirs_overlap_error(fortran_binary): 

    error = get_rirs_overlap_error(fortran_binary) 
    
    ref_error = 0.000869266419

    assert np.allclose(error, ref_error, atol=1.e-6)

