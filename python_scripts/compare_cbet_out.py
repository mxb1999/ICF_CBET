from typing import Iterable, Union
import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import deque
import tables


class Analysis:
    def __init__(self, folder=os.path.abspath('.'), err_analysis='relative'):
        assert os.path.isdir(folder)
        self.directory = folder
        self.err_opt = err_analysis

    def gather_files(path_pattern: str) -> list:
        """Gathers all HDF5 files with the matching prefix

        Args:
            path_pattern (str): the path prefix to be matched

        Returns:
            list: List of open file objects to be handled by cupy operations
        """
        result = deque()
        for f in os.scandir():
            if f.endswith(".h5") and f.beginswith(path_pattern):
                result.append(tables.open_file(f, mode='r'))
        return result

    def compare_files(base_file: tables.File,
                      otherfile: tables.File, param) -> np.ndarray:
        arr1 = np.array(base_file.root[param])
        arr2 = np.array(otherfile.root[param])
        cuarr1 = cp.array(arr1)
        cuarr2 = cp.array(arr2)
        curesult = cp.abs(cp.divide(cp.subtract(cuarr1, cuarr2), cuarr1))
        return cp.asnumpy(curesult)



    def compare(name: Union[str, Iterable]):
        if isinstance(name, Iterable):
            pass
        else:
            pass



