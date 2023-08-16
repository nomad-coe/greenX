"""
! **************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
"""
from pathlib import Path


def find_test_binary(root, binary_name: str) -> Path:
    """ Recursively search the project for a `binary_name`.

    Reference:
    https://stackoverflow.com/questions/18394147/how-to-do-a-recursive-sub-folder-search-and-return-files-in-a-list
    This also has a generator option

    NOTE. This will get slow as one adds many tests.
    Would be better to pass the binary location to the test, however
    one has to deal with pytest's own set of command line args.

    :param root: Root directory in which to begin recursive search.
    :param binary_name: Binary name.
    :return: binary name prepended by its full path.
    """
    found = [path for path in Path(root).rglob(binary_name) if path.is_file()]
    assert len(found) == 1, f"No binary, or multiple binaries were found with the same name: {found}" \
                            f"Began searching from root = {root}"
    return found[0]
