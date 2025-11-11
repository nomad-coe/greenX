"""Testing minimax sine duality error 
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
                          args=['-emin', f'{e_trans.min}', '-emax', f'{e_trans.max}', '-sine'])
    results = runner.run()
    assert results.success, f"Execution of {fortran_binary} failed"
   
    #Parse results
    duality_table = "minimax_duality_error.dat"
    
    return np.genfromtxt(duality_table) 


def test_gx_minimax_duality_error(fortran_binary):

    e_trans = ETrans(0.4, 100.0)

    ref_errors = np.array([[6, 2.50000000E+02, 5.56242714E-02, 2.36021543E-01, 8.74792753E-01],
                           [8, 2.50000000E+02, 1.49383520E-02, 1.11607214E-01, 8.80901314E-01],
                           [10, 2.50000000E+02, 3.87442807E-03, 2.61545954E-02, 9.62458551E-01],
                           [12, 2.50000000E+02, 9.86318209E-04, 1.36267244E-02, 9.39611514E-01],
                           [14, 2.50000000E+02, 2.46521564E-04, 4.19391179E-03, 9.55740373E-01],
                           [16, 2.50000000E+02, 7.70304736E-05, 5.97815737E-04, 9.42329775E-01],
                           [18, 2.50000000E+02, 1.94942501E-05, 1.39775660E-04, 9.45743358E-01],
                           [20, 2.50000000E+02, 4.88850816E-06, 3.73426583E-05, 9.40568100E-01],
                           [22, 2.50000000E+02, 2.68402244E-07, 2.38208860E-05, 1.01253368E+00],
                           [24, 2.50000000E+02, 1.31534818E-07, 2.01595323E-06, 2.04590429E+00],
                           [26, 2.50000000E+02, 3.78405247E-08, 6.37343151E-07, 8.82902975E+00],
                           [28, 2.50000000E+02, 4.30436009E-08, 8.53519608E-07, 1.59367298E+05],
                           [30, 2.50000000E+02, 5.14232643E-08, 1.07148550E-06, 1.05986563E+06],
                           [32, 2.50000000E+02, 4.94132300E-08, 1.02545263E-06, 1.12950579E+06],
                           [34, 2.50000000E+02, 6.66088497E-08, 1.45835462E-06, 5.71640192E+06]])

    tabulated_errors = get_tabulated_errors(fortran_binary, e_trans)

    assert np.allclose(tabulated_errors[:, Column.MaxErrFreqToTime.value],
                       ref_errors[:, Column.MaxErrFreqToTime.value],
                       atol=1.e-7)

    assert np.allclose(tabulated_errors[:, Column.MaxErrTimeToFreq.value],
                       ref_errors[:, Column.MaxErrTimeToFreq.value],
                       atol=1.e-7)

    assert np.allclose(tabulated_errors[0:9, Column.DualityError.value],
                      ref_errors[0:9, Column.DualityError.value],
                      atol=1.e-7)

