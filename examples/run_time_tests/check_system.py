"""
Script for checking system settings.
"""
import os
import platform
import sys
import scipy
import numpy
import astropy
import perf

import MulensModel


print("System: " + os.name + " " + platform.system() + " " + platform.release())
print("Python: " + sys.version)
print("scipy: " + scipy.__version__)
print("numpy: " + numpy.__version__)
print("astropy: " + astropy.__version__)
print("perf: " + perf.__version__)
print("MulensModel: " + MulensModel.__version__)
print("")

# Values below don't matter:
meta = perf.Run([1.0, 1.5, 2.0], warmups=[(1, 3.0)]).get_metadata() 

keys = ['date', 'platform', 'hostname', 'cpu_model_name', 'cpu_count', 
    'python_version']
for key in keys:
    print(key, meta[key])

