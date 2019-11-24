import numpy as np
import os

import MulensModel as mm


SAMPLE_FILE = os.path.join(mm.DATA_PATH, 'unit_test_files', 'FSPL_test_1.dat')


def get_file_params(filename):
    """Read in the model parameters used to create the file"""
    with open(filename) as data_file:
        lines = data_file.readlines()
        ulens_params = lines[2].split()
    return (
        mm.ModelParameters(
            {'t_0': float(ulens_params[1]), 'u_0': float(ulens_params[2]),
             't_E': float(ulens_params[3]), 'rho': float(ulens_params[4])}),
        float(ulens_params[5]))


def get_variables():
    """return a few variables used by 4 test functions below"""
    if 'out' not in get_variables.__dict__:
        names = ['Time', 'b_0', 'b_1', 'Mag_FS', 'Mag_LD', 'Mag']
        data = np.genfromtxt(SAMPLE_FILE, names=names)
        (parameters, gamma) = get_file_params(SAMPLE_FILE)

        point_lens = mm.PointLens(parameters=parameters)
        tau = (data['Time'] - parameters.t_0) / parameters.t_E
        u = np.sqrt(parameters.u_0**2 + tau**2)
        z = u / parameters.rho
        pspl_magnification = (u**2 + 2.) / (u * np.sqrt(u**2 + 4.))

        get_variables.out = (data, gamma, point_lens, u, z, pspl_magnification)
    return get_variables.out


def test_B_0_function():
    """test private _B_0_function"""
    (data, _, point_lens, _, z, _) = get_variables()
    test_b_0 = point_lens._B_0_function(z)
    np.testing.assert_almost_equal(test_b_0, data['b_0'], decimal=5)


def test_B_1_function():
    """test private _B_1_function"""
    (data, _, point_lens, _, z, _) = get_variables()
    test_b_1 = point_lens._B_1_function(z)
    np.testing.assert_almost_equal(test_b_1, data['b_1'], decimal=4)


def test_get_point_lens_finite_source_magnification():
    """test PLFS"""
    (data, _, point_lens, u, _, pspl_magnification) = get_variables()
    test_FSPL = point_lens.get_point_lens_finite_source_magnification(
        u, pspl_magnification)
    np.testing.assert_almost_equal(test_FSPL, data['Mag_FS'], decimal=5)


def test_get_point_lens_limb_darkening_magnification():
    """test PLFS+LD"""
    (data, gamma, point_lens, u, _, pspl_magnification) = get_variables()
    test_FSPL_LD = point_lens.get_point_lens_limb_darkening_magnification(
        u, pspl_magnification, gamma)
    np.testing.assert_almost_equal(test_FSPL_LD/data['Mag_LD'], 1., decimal=4)
