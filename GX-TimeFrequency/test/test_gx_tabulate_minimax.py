import numpy as np
from pathlib import Path
import os
import pytest

from pygreenx.run import BinaryRunner, BuildType, ProcessResults
from pygreenx.utilities import find_test_binary

binary_name = "gx_tabulate_grids.exe"


def parse_std_out(results: ProcessResults):
    """ Parse std.out of gx_tabulate_grids.exe

    Note, if the output format changes in ANY way, this breaks.
    One should switch to structured data.

    :param results subprocess results container.
    :return: Numpy array of errors.
    """
    # Data lines in output
    header = 5
    n_grids = 15
    n_columns = 7

    # Parse std out - only interested in the numerical data
    output_list = results.stdout.decode("utf-8").split('\n')[header:header + n_grids]

    tabulated_errors = np.empty(shape=(n_grids, n_columns))
    for i, line in enumerate(output_list[:]):
        # print(line)
        tabulated_errors[i, :] = np.asarray([float(x) for x in line.split()])

    return tabulated_errors


def test_gx_minimax_grid():

    # NOTE, some of the reference numbers are exceedingly small ~ machine precision
    # and others are very large, which clearly need to be tested with relative tolerance.
    # One could also split up into columns, and test different columns separately.
    ref_errors_small_tau = np.array([[ 6,  0,  2.51200821E-04,  1.88667240E-03,  6.86380679E-03,  5.56242714E-02,  2.50000000E+02],
                                     [ 8,  0,  9.26518238E-04,  7.02310800E-04,  8.97132915E-04,  1.49383520E-02,  2.50000000E+02],
                                     [10,  0,  2.93904398E-02,  4.15097815E-04,  9.85375147E-04,  3.87442807E-03,  2.50000000E+02],
                                     [12,  0,  1.22980779E-03,  5.40966369E-05,  3.23113268E-05,  9.86318209E-04,  2.50000000E+02],
                                     [14,  0,  1.23655331E-03,  1.47661838E-05,  6.49328520E-06,  2.46521564E-04,  2.50000000E+02],
                                     [16,  0,  3.23839478E-03,  4.75444516E-06,  7.88188993E-07,  7.70304736E-05,  2.50000000E+02],
                                     [18,  0,  6.23300175E-03,  1.34768465E-06,  3.45186001E-07,  1.94942501E-05,  2.50000000E+02],
                                     [20,  0,  2.62469243E-03,  2.86296047E-07,  3.67016067E-08,  4.88850816E-06,  2.50000000E+02],
                                     [22,  0,  9.17541725E-02,  1.56304524E-07,  1.11589701E-07,  1.21503959E-06,  2.50000000E+02]])

    ref_errors_big_tau = np.array([[24,  0,  5.95421450E+01,  9.11762981E-08,  9.65139449E-09,  1.43312081E-06,  2.50000000E+02],
                                   [26,  0,  9.37145464E+02,  3.33178951E-08,  3.39771986E-09,  5.16593492E-07,  2.50000000E+02],
                                   [28,  0,  7.58978121E+06,  2.91854558E-08,  4.27640628E-09,  4.63466294E-07,  2.50000000E+02],
                                   [30,  0,  2.88001056E+08,  3.72208007E-08,  2.57757957E-09,  7.25625128E-07,  2.50000000E+02],
                                   [32,  0,  1.50874511E+08,  2.89721647E-08,  2.68576848E-09,  3.26754726E-07,  2.50000000E+02],
                                   [34,  0,  1.03478980E+08,  3.74825302E-08,  2.22192890E-09,  2.66427793E-07,  2.50000000E+02]])

    # Find binary
    root = Path(__file__).parent.parent.parent
    binary = find_test_binary(root, binary_name)
    assert os.path.isfile(binary), f'Binary {binary_name} cannot be found'

    # Run binary with specified emin and emax
    emin, emax = 0.4, 100.0
    runner = BinaryRunner(binary, BuildType.serial,
                          args=['table', '-emin', f'{emin}', '-emax', f'{emax}'])
    results = runner.run()
    assert results.success, f"Execution of {binary} failed"

    # Compare parsed results to reference values
    tabulated_errors = parse_std_out(results)
    print(tabulated_errors[0:9, :] - ref_errors_small_tau)
    print(tabulated_errors[9:, :] - ref_errors_big_tau)

    assert np.allclose(tabulated_errors[0:9, :], ref_errors_small_tau, atol=1.e-7), \
        'Difference in results with small tau values'

    # TODO Investigate if these larger grids are stable
    # Alex gets different numbers for grid 28 tau between his Mac with openblas,
    # and the Ubuntu CI with blas/lapack
    # assert np.allclose(tabulated_errors[9:, :], ref_errors_big_tau, atol=0.01), \
    #     'Difference in results with large tau values'
