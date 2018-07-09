import numpy as np
from numpy_functions import numpy_chi2_v3


MAG_ZEROPOINT = 18.

file_1 = "test_100.txt"
data_1 = np.loadtxt(file_1, unpack=True)
time = data_1[0]
obs_flux = np.power(10., -0.4 * (data_1[1] - MAG_ZEROPOINT))
obs_flux_err = data_1[2] * obs_flux * 0.4 * np.log(10.)

t_0 = 2456900.
u_0 = 0.01
t_E = 20.
