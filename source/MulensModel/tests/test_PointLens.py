import numpy as np
import os

import MulensModel
from MulensModel import PointLens, ModelParameters

DATA_PATH = os.path.join(MulensModel.MODULE_PATH, 'data')
SAMPLE_FILE = os.path.join(DATA_PATH, 'FSPL_test_1.dat')

def get_file_params(filename):
    """Read in the model parameters used to create the file"""
    with open(filename) as data_file:
        lines = data_file.readlines()
        ulens_params = lines[2].split()
    return (
        ModelParameters(
            {'t_0': float(ulens_params[1]), 'u_0': float(ulens_params[2]), 
             't_E': float(ulens_params[3]), 'rho': float(ulens_params[4])}),
        float(ulens_params[5]))

DATA = np.genfromtxt(
    SAMPLE_FILE, dtype=None, 
    names=['Time', 'b_0', 'b_1', 'Mag_FS', 'Mag_LD', 'Mag'])
(PARAMETERS, GAMMA) = get_file_params(SAMPLE_FILE)

point_lens = PointLens(parameters=PARAMETERS)
tau = (DATA['Time'] - PARAMETERS.t_0) / PARAMETERS.t_E
u = np.sqrt(PARAMETERS.u_0**2 + tau**2)
z = u / PARAMETERS.rho
pspl_magnification = (u**2 + 2.) / (u * np.sqrt(u**2 + 4.))

mid_index = 0

def test_B_0_function():
    test_b_0 = point_lens._B_0_function(z)
    np.testing.assert_almost_equal(test_b_0, DATA['b_0'], decimal=5)

def test_B_1_function():
    print(PARAMETERS)
    test_b_1 = point_lens._B_1_function(z)
    print(
        test_b_1[mid_index], 
        DATA['b_1'][mid_index])
    np.testing.assert_almost_equal(test_b_1, DATA['b_1'], decimal=5)

def test_get_point_lens_finite_source_magnification():
    test_FSPL = point_lens.get_point_lens_finite_source_magnification(
        u, pspl_magnification)
    print(
        test_FSPL[mid_index], 
        DATA['Mag_FS'][mid_index])
    np.testing.assert_almost_equal(test_FSPL, DATA['Mag_FS'], decimal=5)

def test_get_point_lens_limb_darkening_magnification():
    print('gamma = ', GAMMA)
    test_FSPL_LD = point_lens.get_point_lens_limb_darkening_magnification(
        u, pspl_magnification, GAMMA)
    print(
        test_FSPL_LD[mid_index], 
        DATA['Mag_LD'][mid_index])
    np.testing.assert_almost_equal(test_FSPL_LD, DATA['Mag_LD'], decimal=5)


