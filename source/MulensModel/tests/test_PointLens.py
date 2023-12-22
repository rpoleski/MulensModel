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
        trajectory = mm.Trajectory(data['Time'], parameters)

        get_variables.out = (data, gamma, trajectory)
    return get_variables.out


def test_B_0_function():
    """test private _B_0_function"""
    (data, _, trajectory) = get_variables()
    point_lens = mm.FiniteSourceUniformGould94Magnification(
        trajectory=trajectory)
    test_b_0 = point_lens._B_0_function()
    np.testing.assert_almost_equal(test_b_0, data['b_0'], decimal=5)


def test_B_1_function():
    """test private _B_1_function"""
    (data, gamma, trajetory) = get_variables()
    test_FSPL_LD = mm.FiniteSourceLDYoo04Magnification(
        trajectory=trajectory, gamma=gamma)
    test_b_1 =  test_FSPL_LD._B_1_function()
    np.testing.assert_almost_equal(test_b_1, data['b_1'], decimal=4)


def test_get_point_lens_finite_source_magnification():
    """test PLFS"""
    (data, _, trajectory) = get_variables()
    test_FSPL = mm.FiniteSourceUniformGould94Magnification(
        trajectory=trajectory)
    fspl_magnification = test_FSPL.get_magnification()
    np.testing.assert_almost_equal(fspl_magnificationL, data['Mag_FS'], decimal=5)


def test_get_point_lens_limb_darkening_magnification():
    """test PLFS+LD"""
    (data, gamma, trajetory) = get_variables()
    test_FSPL_LD = mm.FiniteSourceLDYoo04Magnification(
        trajectory=trajectory, gamma=gamma)
    fspl_magnification = test_FSPL_LD.get_magnification()
    np.testing.assert_almost_equal(fspl_magnification/data['Mag_LD'], 1., decimal=4)
