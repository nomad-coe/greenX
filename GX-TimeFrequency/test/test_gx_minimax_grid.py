""" Inputs and assertions for test_gx_minimax_grid.f90

Would be better if one had python API to GX, and avoid the need for this intermediate binary.
Indeed, MUCH of the complication regarding running tests comes from this.

Issues and solutions:
---------------------
Problem
* python needs to know where the test binary is located.
Solution
* The robust but less-than-ideal approach solution is to have python recursively
  search through the directories for it.
  This work very well if one defines the tests to be run from a given build
  directory.

Problem
* Where should the python test write the input data for the fortran binary?
Solution
* This is solved using tmp_path, which gives an absolute path.


* Path for passing to the fortran binary     () This is solved using tmp_path, which gives an absolute path
* where to read the result from              () This is solved using tmp_path, which gives an absolute path
"""
import numpy as np
from pathlib import Path

from pygreenx.run import BinaryRunner, BuildType
from pygreenx.utilities import find_test_binary

binary_name = "test_minimax.exe"


def mock_file(tmp_path, inputs_str):
    """ Write an inputs string to file.
    :param tmp_path:
    :param inputs_str:
    :return:
    """
    # Write test inputs to file
    file: Path = tmp_path / "inputs.dat"
    file.write_text(inputs_str)
    assert file.exists(), f"{file.name} does not exist"
    return file


def set_input(tmp_path, n_mesh_points, e_trans_min, e_trans_max):
    """ Set inputs in a string.

    :param tmp_path:
    :param n_mesh_points:
    :param e_trans_min:
    :param e_trans_max:
    :return:
    """
    inputs_template = """
{}
n_mesh_points     {}
e_transition_min  {}
e_transition_max  {}
"""
    inputs_str = inputs_template.format(tmp_path, n_mesh_points, e_trans_min, e_trans_max)
    return inputs_str


def test_gx_minimax_grid(tmp_path):
    # Test inputs
    n_mesh_points = 10
    e_trans_min, e_trans_max = 2., 30.

    # Fortran input file
    inputs_str = set_input(tmp_path, n_mesh_points, e_trans_min, e_trans_max)
    file = mock_file(tmp_path, inputs_str)

    # Run test
    # root MUST be defined in the test script
    # However, not robust to changes in the nesting of directories.
    # For example, if I change <BUILD_DIR>/test to <BUILD_DIR>/test/time-frequency
    # this command requires updating (and vice versa)
    root = Path(__file__).parent.parent.parent
    binary = find_test_binary(root, binary_name)
    runner = BinaryRunner(binary, BuildType.serial, args=[file.as_posix()])
    results = runner.run()
    assert results.success, f"Execution of {binary} failed"

    # Parse results
    tau_file = tmp_path / "tau.dat"
    freq_file = tmp_path / "freq.dat"

    tau_grid_weights = np.genfromtxt(tau_file)
    freq_grid_weights = np.genfromtxt(freq_file)

    ref_tau_grid_weights = np.array([[0.005918, 0.01525399],
                                     [0.03177845, 0.03688566],
                                     [0.08092057, 0.0622535],
                                     [0.15862351, 0.09469842],
                                     [0.27439677, 0.13947465],
                                     [0.44433573, 0.20467764],
                                     [0.69464678, 0.30283334],
                                     [1.06797564, 0.45564411],
                                     [1.63956744, 0.7122558],
                                     [2.58159534, 1.25940329]])

    ref_freq_grid_weights = np.array([[0.37397426, 3.02916016],
                                      [1.17910655, 3.49425362],
                                      [2.16745299, 4.51933688],
                                      [3.50267933, 6.31709192],
                                      [5.42267468, 9.28613879],
                                      [8.30517692, 14.19510809],
                                      [12.81707194, 22.74856946],
                                      [20.35677335, 39.7693943],
                                      [34.86933737, 84.9612311],
                                      [75.79619716, 320.31845324]])

    assert np.allclose(tau_grid_weights, ref_tau_grid_weights)
    assert np.allclose(freq_grid_weights, ref_freq_grid_weights)
