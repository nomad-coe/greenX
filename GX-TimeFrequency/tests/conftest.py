"""
Module to define metafunc entries from command-line arguments.

metafunc can be passed to pytest fixtures.
"""

def pytest_addoption(parser):
    """
    Custom command-line options parsed added to pytest.
    Function used by pytest.

    These are ONLY used in the test runner to determine the full path to
    the fortran executable and the run settings.
    """

    parser.addoption('--root',
                     type=str,
                     dest='root',
                     required=True,
                     help='Full path for GreenX root'
                     )

    parser.addoption('--binary',
                     type=str,
                     dest='binary',
                     required=True,
                     help='Name of the executable with full path prepending'
                     )


def pytest_generate_tests(metafunc):
    """
    Add parameter or parameters to metafunc.
    Function used by pytest.

    :param metafunc: Look me up
    """
    for option in ['root', 'binary']:
        if option in metafunc.fixturenames:
            metafunc.parametrize(option, metafunc.config.getoption(option))