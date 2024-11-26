import os
from glob import glob
import ctypes

_library_file = glob(os.path.dirname(__file__) + "/libcgasdynamicslib.*")[0]
_library = ctypes.CDLL(_library_file)

_library.gs_test.argtypes = ()
