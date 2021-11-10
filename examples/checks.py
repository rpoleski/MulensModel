import os
import platform
print("System: " + os.name + " " + platform.system() + " " + platform.release())

import sys
print("Python: " + sys.version)

import ctypes
print("ctypes: " + ctypes.__version__)

import scipy
print("scipy: " + scipy.__version__)

import numpy
print("numpy: " + numpy.__version__)

import matplotlib
print("matplotlib: " + matplotlib.__version__)

import astropy
print("astropy: " + astropy.__version__)

import sympy
print("sympy: " + sympy.__version__)

import MulensModel
print("MulensModel: " + MulensModel.__version__)

