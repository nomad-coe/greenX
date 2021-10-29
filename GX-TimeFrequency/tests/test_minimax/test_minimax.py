"""
Minimax frequency grid regression tests

Run from GreenX's root:

pytest -s GX-TimeFrequency/tests/test_minimax/test_minimax.py --root <GREENX_ROOT> --binary path/to/<MINIMAX_TEST.EXE>

"""
import pytest
import numpy as np
from pathlib import Path

#from pygreenx.run import BinaryRunner, BuildType


@pytest.fixture(scope="session")
def paths(pytestconfig):
    greenx_root = Path(pytestconfig.getoption("root"))
    binary = Path(pytestconfig.getoption('binary'))
    test_dir = greenx_root / f'GX-TimeFrequency/tests/test_minimax'
    return binary, test_dir


# TODO(Alex) Generate reference data for each grid size
@pytest.mark.parametrize("grid_size", ['32'])
def test_minimax_grids_h20_molecule(paths,  grid_size):
    """
    Test all minimax grids for a eigenvalues of an H2O molecule 
    """
    binary, test_dir = paths

    # Input file
    eigenvalue_input = test_dir / f'inputs/H2O_eigenvalues.txt'
    homo_index = '5'

    # Output grid
    ref_file = test_dir / f'refs/{grid_size}_freq_points.dat'
    output_file = test_dir / f'output_{grid_size}_freq_points.dat'

    # Run fortran code
    minimax_args = ['-s', grid_size, '-e', eigenvalue_input, '-h', homo_index, '-output', output_file]
    test_case = BinaryRunner(binary, BuildType.serial, args=minimax_args)
    results = test_case.run()
    assert results.success, "minimax test binary execution failed"

    # Compare output to reference 
    ref = np.loadtxt(ref_file, skiprows=2)
    output = np.loadtxt(output_file, skiprows=2)
    assert np.allclose(ref, output, atol=1.e-8), 'Reference minimax grid disagrees with output'
