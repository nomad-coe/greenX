from pathlib import Path
import pytest
import os
from typing import Union


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
                     required=False,
                     help='Full path for GreenX root'
                     )

    parser.addoption('--binary',
                     type=str,
                     dest='binary',
                     required=False,
                     help='Name of the executable with full path prepending'
                     )


# Define metafunc entries from command-line arguments.
# metafunc can be passed to pytest fixtures.
def pytest_generate_tests(metafunc):
    """ Add parameter or parameters to metafunc.

    :param metafunc
    """
    for option in ['root', 'binary']:
        if option in metafunc.fixturenames:
            metafunc.parametrize(option, metafunc.config.getoption(option))


@pytest.fixture(scope="session")
def greenx_build_root() -> str:
    """ Get GreenX build directory from an environment variable
    set by CMake, during the build process.

    :return: Environment variable string
    """
    ENV_VAR = 'GX_BUILD_DIR'
    return os.getenv(ENV_VAR)


@pytest.fixture(scope="session")
def all_binaries(greenx_build_root, valid_binary_extensions=None) -> dict:
    """ List all binaries recursively found in a directory.

    Given some top-level directory, recursively find all binaries present,
    where binary is defined according to valid_binary_extensions.

    :param greenx_build_root: Build directory for GreenX
    :param valid_binary_extensions:
    :return: result with {key:value} = {file_name: Path/to/file}
    """
    assert greenx_build_root is not None, 'GX build directory cannot be found. ' \
                                          'Try `export GX_BUILD_DIR=<PATH/2/BUILD>`'

    if valid_binary_extensions is None:
        valid_binary_extensions = ['exe', 'x']

    binaries = []
    for ext in valid_binary_extensions:
        matches = Path(greenx_build_root).rglob('*' + ext)
        binaries += [path for path in matches if path.is_file()]

    # Pack for convenient look-up
    result = {}
    for binary in binaries:
        binary: Path
        result[binary.name] = binary.parent

    return result


@pytest.fixture()
def get_binary(all_binaries):
    def wrapper(name: str) -> Union[Path, None]:
        """ Helper function to retrieve absolute binary path.

        :param name: Binary name (including extension)
        :return: Binary prepended by absolute path, or None if the binary
        is not present in all_binaries.
        """
        path: Path = all_binaries.get(name)

        if path is None:
            return path

        return path / name
    return wrapper

