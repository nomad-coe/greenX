from pathlib import Path
import pytest
import os
from typing import Union


def pytest_addoption(parser):
    """ Custom command-line options parsed added to pytest.
    """
    parser.addoption('--root',
                     type=str,
                     dest='root',
                     required=False,
                     help='Full path for GreenX root'
                     )

    parser.addoption('--build-dir',
                     type=str,
                     dest='build_dir',
                     default=None,
                     required=False,
                     help='Full path to build directory'
                     )

    parser.addoption('--binary',
                     type=str,
                     dest='binary',
                     required=False,
                     help='Name of the executable with full path prepended'
                     )


@pytest.fixture(scope="session")
def greenx_build_root(request) -> str:
    """ Get GreenX build directory.

    :return: Environment variable string
    """
    # First, try from command line arg.
    build_dir = request.config.getoption('build_dir')
    if build_dir:
        if not Path(build_dir).is_dir():
            raise NotADirectoryError(f'{build_dir} is not a directory')
        return build_dir

    #  If this is not passed, look for an environment variable
    if build_dir is None:
        return os.getenv('GX_BUILD_DIR')


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


@pytest.fixture()
def get_binary_from_build(get_binary, greenx_build_root):
    """ Get a binary prepended by absolute path, from the GX build directory.

    :param get_binary: pytest fixture to retrieve absolute binary path.
    :param greenx_build_root: GX root directory.
    :return:
    """
    def wrapper(name: str):
        _binary = get_binary(name)
        assert _binary is not None, f'{name} cannot be found in path {greenx_build_root}'
        print(f'Binary source: {_binary}')
        return _binary
    return wrapper


# Run automatically
@pytest.fixture(autouse=True)
def testdir(tmp_path: Path):
    """ Test directory fixture.

    Provide a fixture to switch to tmp_path as a working test directory.

    :param tmp_path: Temporary directory unique to the test invocation
    """
    oldpath = os.getcwd()
    os.chdir(tmp_path.as_posix())
    yield
    os.chdir(oldpath)
    print(f'Testdir: {tmp_path.as_posix()}')
