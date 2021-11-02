"""
Module to run an executable
"""
import os
import pathlib
import enum
import subprocess
from typing import Optional, Union, List


@enum.unique
class BuildType(enum.Enum):
    serial = enum.auto()
    pure_mpi = enum.auto()
    pure_omp = enum.auto()
    omp_mpi = enum.auto()


default_run_cmd = {BuildType.serial: './',
                   BuildType.pure_mpi: 'mpirun',
                   BuildType.pure_omp: './',
                   BuildType.omp_mpi: 'mpirun'
                   }


class ProcessResults:
    def __init__(self, stdout, stderr, return_code: int):
        self.stdout = stdout
        self.stderr = stderr
        self.return_code = return_code
        self.success = return_code == 0


class BinaryRunner:
    def __init__(self,
                 binary: str,
                 build_type: BuildType,
                 omp_num_threads: Optional[int] = 1,
                 mpi_processes: Optional[int] = 1,
                 run_cmd: Optional[Union[str, dict]] = default_run_cmd,
                 args: Optional[List[str]] = None
                 ) -> None:

        self.binary = binary
        self.build_type = build_type
        self.omp_num_threads = omp_num_threads
        self.mpi_processes = mpi_processes
        self.run_cmd = run_cmd[self.build_type] if isinstance(run_cmd, dict) \
            else run_cmd
        self.args = [''] if args is None else args

        try:
            pathlib.Path(self.binary).resolve(strict=True)
        except FileNotFoundError:
            raise FileNotFoundError("Binary not found:" + self.binary)

        assert omp_num_threads > 0, "Number of OMP threads must be > 0"

        assert mpi_processes > 0, "Number of MPI processes must be > 0"

    def get_run_list(self) -> list:
        """
        Generate a list of run arguments for subprocess.run

        TODO(Alex) Extend to pass MPI arguments
        """
        if self.build_type in [BuildType.serial, BuildType.pure_omp]:
            return [self.binary] + self.args

        elif self.build_type in [BuildType.pure_mpi, BuildType.omp_mpi]:
            return [self.run_cmd, '-np', self.mpi_processes,
                    self.binary] + self.args

        else:
            raise ValueError('Build type not recognised')

    def run(self, run_args: Optional[list] = None) -> ProcessResults:
        """
        Run a binary

        :param Optional[list] run_args: List of arguments required by subprocess.run]. 
                                        Defaults to None.
        """
        if run_args is None:
            run_args = self.get_run_list()

        my_env = {**os.environ, "OMP_NUM_THREADS": str(self.omp_num_threads)}
        result = subprocess.run(run_args, env=my_env, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        return ProcessResults(result.stdout, result.stderr, result.returncode)
