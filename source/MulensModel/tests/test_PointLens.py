import numpy as np

import MulensModel.MODULE_PATH
from MulensModel import PointLens

DATA_PATH = os.path.join(MulensModel.MODULE_PATH, 'data')
SAMPLE_FILE = os.path.join(DATA_PATH, 'FSPL_test.dat')

def get_file_params(filename):
    """Read in the model parameters used to create the file"""
    with open(filename) as data_file:
        lines = data_file.readlines()
        ulens_params = lines[3].split()
    return (
        {'t_0': ulens_params[1], 'u_0': ulens_params[2], 
         't_E': ulens_params[3], 'rho': ulens_params[3]},
        gamma)

DATA = np.genfromtxt(
    SAMPLE_FILE, dtype=None, 
    names=['Time', 'b0', 'b1', 'Mag_FS', 'Mag_LD', 'Mag'])

(PARAMETERS, GAMMA) = get_file_params(SAMPLE_FILE)

