"""Testing minimax sine duality error with tikhonov regularization  
   Copyright (C) 2020-2025 GreenX library
   This file is distributed under the terms of the APACHE2 License.
"""
import enum
import numpy as np
import pytest

from pygreenx.run import BinaryRunner, BuildType, ProcessResults


@pytest.fixture()
@pytest.mark.usefixtures("get_binary", "greenx_build_root")
def fortran_binary(get_binary, greenx_build_root):
    name = 'test_gx_minimax_duality_error.exe'
    _binary = get_binary(name)
    assert _binary is not None, f'{name} cannot be found in {greenx_build_root}'
    print(f'Binary source: {_binary}')
    return _binary


class ETrans:
    def __init__(self, min, max):
        self.min = min
        self.max = max

class Column(enum.Enum):
    """Reference Data Column Labels """
    NumPoints = 0
    ERatio = 1
    MaxErrFreqToTime = 2
    MaxErrTimeToFreq = 3
    DualityError = 4


def get_tabulated_errors(fortran_binary, e_trans: ETrans):
    runner = BinaryRunner(fortran_binary, BuildType.serial,
                          args=['-emin', f'{e_trans.min}', '-emax', f'{e_trans.max}', '-sine', '-regularization', '0.02'])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"
   
    #Parse results
    duality_table = "minimax_duality_error.dat"
    
    return np.genfromtxt(duality_table) 


def test_gx_minimax_duality_error(fortran_binary):

    e_trans = ETrans(0.4, 100.0)

    ref_errors = np.array([[6, 2.50000000E+02, 6.37542287E-02, 2.68482277E-01, 8.47611913E-01],
                           [8, 2.50000000E+02, 3.34884766E-02, 2.14747858E-01, 9.45288283E-01],
                           [10, 2.50000000E+02, 3.07259007E-02, 1.73419958E-01, 1.04509614E+00],
                           [12, 2.50000000E+02, 2.95633756E-02, 2.07312027E-01, 1.12233487E+00],
                           [14, 2.50000000E+02, 2.78894111E-02, 2.17166449E-01, 1.05395626E+00],
                           [16, 2.50000000E+02, 2.66005624E-02, 2.40660950E-01, 1.01303721E+00],
                           [18, 2.50000000E+02, 2.52503847E-02, 2.51167867E-01, 1.00651221E+00],
                           [20, 2.50000000E+02, 2.40563752E-02, 2.63992682E-01, 1.00018000E+00],
                           [22, 2.50000000E+02, 2.26287888E-02, 2.46955881E-01, 9.92094591E-01],
                           [24, 2.50000000E+02, 2.25973528E-02, 2.91937768E-01, 9.90679995E-01],
                           [26, 2.50000000E+02, 2.23747422E-02, 2.98938299E-01, 9.92758897E-01],
                           [28, 2.50000000E+02, 2.25060973E-02, 2.97170880E-01, 9.92244263E-01],
                           [30, 2.50000000E+02, 2.25094063E-02, 2.95860549E-01, 9.97887535E-01],
                           [32, 2.50000000E+02, 2.25550485E-02, 2.96175702E-01, 9.99375316E-01],
                           [34, 2.50000000E+02, 2.26898031E-02, 2.94091988E-01, 9.99869417E-01]])

    tabulated_errors = get_tabulated_errors(fortran_binary, e_trans)

    assert np.allclose(tabulated_errors[:, Column.MaxErrFreqToTime.value],
                       ref_errors[:, Column.MaxErrFreqToTime.value],
                       atol=1.e-7)

    assert np.allclose(tabulated_errors[:, Column.MaxErrTimeToFreq.value],
                       ref_errors[:, Column.MaxErrTimeToFreq.value],
                       atol=1.e-7)

    assert np.allclose(tabulated_errors[:, Column.DualityError.value],
                      ref_errors[:, Column.DualityError.value],
                      atol=1.e-7)

