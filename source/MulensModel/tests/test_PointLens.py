import numpy as np
import os


import MulensModel as mm
from test_FitData import create_0939_parallax_model, SAMPLE_FILE_03, \
    SAMPLE_FILE_03_EPH

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
    (data, gamma, trajectory) = get_variables()
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
    np.testing.assert_almost_equal(
        fspl_magnification, data['Mag_FS'], decimal=5)


def test_get_point_lens_limb_darkening_magnification():
    """test PLFS+LD"""
    (data, gamma, trajectory) = get_variables()
    test_FSPL_LD = mm.FiniteSourceLDYoo04Magnification(
        trajectory=trajectory, gamma=gamma)
    fspl_magnification = test_FSPL_LD.get_magnification()
    np.testing.assert_almost_equal(
        fspl_magnification/data['Mag_LD'], 1., decimal=4)

def test_fspl_noLD():
    """
    check if FSPL magnification is calculate properly
    """
    t_0 = 2456789.012345
    t_E = 23.4567
    u_0 = 1e-4
    rho = 1e-3
    t_vec = np.array([-(rho**2-u_0**2)**0.5, 0., ((0.5*rho)**2-u_0**2)**0.5])
    t_vec = t_vec * t_E + t_0

    params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})

    trajectory = mm.Trajectory(t_vec, params)

    mag_curve = mm.FiniteSourceUniformGould94Magnification(trajectory=trajectory)
    results = mag_curve.get_magnification()

    u = np.array([rho, u_0, 0.5*rho])
    pspl = (u**2 + 2.) / np.sqrt(u**2 * (u**2 + 4.))
    expected = np.array([1.27323965, 0.19949906, 0.93421546])
    # These values were calculated by Andy Gould (file b0b1.dat).
    expected *= pspl

    np.testing.assert_almost_equal(expected, results, decimal=4)

def test_get_d_u_d_params():
    """
    Test that calculating derivatives with an ephemeris file is different from
    without an ephemeris file.
    """
    parameters = ['pi_E_N', 'pi_E_E']
    model_with_par = create_0939_parallax_model()
    data_ephm = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    parallax = {'earth_orbital': True,
                'satellite': True,
                'topocentric': True}

    traj_ephm = mm.Trajectory(
        data_ephm.time, parameters=model_with_par.parameters,
        satellite_skycoord=data_ephm.satellite_skycoord,
        coords=model_with_par.coords, parallax=parallax)
    pl_ephm = mm.PointSourcePointLensMagnification(traj_ephm)
    derivs_ephm = pl_ephm.get_d_u_d_params(parameters)

    traj_no_ephm = mm.Trajectory(
        data_ephm.time, parameters=model_with_par.parameters,
        coords=model_with_par.coords, parallax=parallax)
    pl_no_ephm = mm.PointSourcePointLensMagnification(traj_no_ephm)
    derivs_no_ephm = pl_no_ephm.get_d_u_d_params(parameters)

    for param in parameters:
        ratio = derivs_ephm[param] / derivs_no_ephm[param]
        assert (np.abs(ratio - 1.) > 0.001).all()
