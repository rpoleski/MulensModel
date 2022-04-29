import unittest
import numpy as np

from astropy import units as u

import MulensModel as mm


class TestModelParameters(unittest.TestCase):
    def test_too_many_parameters_for_init(self):
        with self.assertRaises(KeyError):
            mp = mm.ModelParameters(
                {'pi_E': (1., 1.), 'pi_E_N': 1.})
        with self.assertRaises(KeyError):
            mp = mm.ModelParameters(
                {'pi_E': (1., 1.), 'pi_E_E': 1.})


def test_init_parameters():
    """
    make sure that parameters are properly stored
    """
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E.value)


def test_repr_parameters():
    """
    Test str(ModelParameters) or __repr__(ModelParameters)
    """
    t_0 = 2456141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

    out_1 = "    t_0 (HJD)       u_0    t_E (d) \n"
    out_2 = "2456141.59300  0.542500    62.6300 \n"

    assert (out_1 + out_2) == str(params)


def test_rho_t_e_t_star():
    """check if conversions between rho, t_E, and t_star work ok"""
    t_0 = 2450000.
    u_0 = 0.1
    t_E_1 = 20.
    t_E_2 = t_E_1 * u.day
    rho = 0.001
    t_star_1 = t_E_1 * rho
    t_star_2 = t_E_2 * rho

    for (t_E, t_star) in zip([t_E_1, t_E_2], [t_star_1, t_star_2]):
        params_1 = mm.ModelParameters(
            {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})
        np.testing.assert_almost_equal(params_1.t_star, t_star_1)

        params_2 = mm.ModelParameters(
            {'t_0': t_0, 'u_0': u_0, 't_star': t_star, 'rho': rho})
        np.testing.assert_almost_equal(params_2.t_E, t_E_1)

        params_3 = mm.ModelParameters(
            {'t_0': t_0, 'u_0': u_0, 't_star': t_star, 't_E': t_E})
        np.testing.assert_almost_equal(params_3.rho, rho)


class test(unittest.TestCase):
    def test_too_much_rho_t_e_t_star(self):
        with self.assertRaises(KeyError):
            t_0 = 2450000.
            u_0 = 0.1
            t_E = 20. * u.day
            rho = 0.001
            t_star = t_E * rho
            mm.ModelParameters({
                't_0': t_0, 'u_0': u_0, 't_E': t_E,
                'rho': rho, 't_star': t_star})


def test_update():
    """
    check if the number of parameters can be changed
    """
    t_0 = 2456141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    params.as_dict().update({'rho': 0.001})

    assert len(params.parameters.keys()) == 4


def test_orbital_motion_1():
    """basic tests of orbital motion"""
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 30.}
    dict_motion = dict_static.copy()
    dict_motion.update({'dalpha_dt': 10., 'ds_dt': 0.1})
    static = mm.ModelParameters(dict_static)
    motion = mm.ModelParameters(dict_motion)

    assert static.is_static()
    assert not motion.is_static()

    epoch_1 = dict_static['t_0'] - 18.2625
    epoch_2 = dict_static['t_0']
    epoch_3 = dict_static['t_0'] + 18.2625

    # Test get_s() and get_alpha() for static case.
    assert static.get_s(epoch_1) == static.get_s(epoch_2)
    assert static.get_s(epoch_1) == static.get_s(epoch_3)
    assert static.get_s(epoch_1) == dict_static['s']
    assert static.get_alpha(epoch_1) == static.get_alpha(epoch_3)
    assert static.get_alpha(epoch_1) == dict_static['alpha'] * u.deg

    # Test get_s() and get_alpha() for orbital motion case.
    np.testing.assert_almost_equal(motion.get_alpha(epoch_1).value, 29.5)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_2).value, 30.)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_3).value, 30.5)
    np.testing.assert_almost_equal(motion.get_s(epoch_1), 1.2295)
    np.testing.assert_almost_equal(motion.get_s(epoch_2), 1.2345)
    np.testing.assert_almost_equal(motion.get_s(epoch_3), 1.2395)

    # Test arguments as list or array.
    np.testing.assert_almost_equal(
        motion.get_alpha([epoch_1, epoch_2, epoch_3]).value,
        [29.5, 30., 30.5])
    np.testing.assert_almost_equal(
        motion.get_alpha(np.array([epoch_1, epoch_2, epoch_3])).value,
        [29.5, 30., 30.5])
    np.testing.assert_almost_equal(
        motion.get_s([epoch_1, epoch_2, epoch_3]),
        [1.2295, 1.2345, 1.2395])
    np.testing.assert_almost_equal(
        motion.get_s(np.array([epoch_1, epoch_2, epoch_3])),
        [1.2295, 1.2345, 1.2395])

    # Test get_alpha() units.
    assert static.get_alpha(epoch_1).unit == u.deg
    assert static.get_alpha(epoch_2).unit == u.deg
    assert motion.get_alpha(epoch_1).unit == u.deg
    assert motion.get_alpha(epoch_2).unit == u.deg

    assert motion.alpha == 30. * u.deg
    assert motion.s == 1.2345


def test_t_0_kep():
    """test reference epoch for orbital motion"""
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 30.}
    dict_motion = dict_static.copy()
    dict_motion.update({'dalpha_dt': 10., 'ds_dt': 0.1})
    dict_motion_2 = dict_motion.copy()
    dict_motion_2.update({'t_0_kep': 2456789.01234-18.2625})
    motion = mm.ModelParameters(dict_motion)
    motion_2 = mm.ModelParameters(dict_motion_2)

    assert motion.t_0_kep == motion.t_0
    assert motion_2.t_0_kep == dict_motion_2['t_0_kep']

    epoch_1 = dict_static['t_0'] - 18.2625
    epoch_2 = dict_static['t_0']

    # Test motion.
    np.testing.assert_almost_equal(motion.get_alpha(epoch_1).value, 29.5)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_2).value, 30.)
    np.testing.assert_almost_equal(motion.get_s(epoch_1), 1.2295)
    np.testing.assert_almost_equal(motion.get_s(epoch_2), 1.2345)

    # Test motion_2.
    np.testing.assert_almost_equal(motion_2.get_alpha(epoch_1).value, 30.)
    np.testing.assert_almost_equal(motion_2.get_alpha(epoch_2).value, 30.5)
    np.testing.assert_almost_equal(motion_2.get_s(epoch_1), 1.2345)
    np.testing.assert_almost_equal(motion_2.get_s(epoch_2), 1.2395)


def test_orbital_motion_gammas():
    """test .gamma_parallel .gamma_perp .gamma"""
    dict_params = {'t_0': 2457123.456, 'u_0': 0.0345, 't_E': 30.00,
                   's': 1.5, 'ds_dt': 0.5, 'q': 0.987,
                   'alpha': 12.345, 'dalpha_dt': 50.}
    params = mm.ModelParameters(dict_params)

    # Test values.
    np.testing.assert_almost_equal(params.gamma_parallel.value, 0.333333333)
    np.testing.assert_almost_equal(params.gamma_perp.value, -0.872664626)
    np.testing.assert_almost_equal(params.gamma.value, 0.934159869)

    # Test units.
    assert params.gamma_parallel.unit == 1. / u.year
    assert params.gamma_perp.unit == u.rad / u.year
    assert params.gamma.unit == 1. / u.year


def test_binary_source():
    """
    Test if binary source parameters are properly set and changed
    """
    params = mm.ModelParameters({
        't_0_1': 2455000.1, 'u_0_1': 0.05, 't_star_1': 0.025,
        't_0_2': 2455123.4, 'u_0_2': 0.15, 't_star_2': 0.050,
        't_E': 25.})
    assert params.t_E == 25.
    assert params.source_1_parameters.rho == 0.001
    assert params.source_2_parameters.rho == 0.002

    params.t_0_1 = 2456789.012345
    assert params.t_0_1 == 2456789.012345
    assert params.source_1_parameters.t_0 == 2456789.012345

    params.t_star_2 = 0.075
    assert params.source_2_parameters.t_star == 0.075
    assert params.t_star_2 == 0.075
    assert params.source_2_parameters.rho == 0.003
    assert params.rho_2 == 0.003

    params.t_E = 12.5
    assert params.t_star_1 == 0.025
    assert params.source_1_parameters.t_star == 0.025
    assert params.rho_1 == 0.002
    assert params.source_1_parameters.rho == 0.002
    assert params.t_star_2 == 0.075
    assert params.source_2_parameters.t_star == 0.075
    assert params.rho_2 == 0.006
    assert params.source_2_parameters.rho == 0.006


def test_n_lenses():
    """
    Test if lenses are counted properly
    """
    p_1 = mm.ModelParameters({'t_0': 10, 'u_0': 1, 't_E': 3, 'rho': 0.001})
    p_2 = mm.ModelParameters({'t_0': 10, 'u_0': 1, 't_E': 3, 'rho': 0.001,
                              's': 10, 'q': 0.5, 'alpha': 100.})
    p_3 = mm.ModelParameters({'x_caustic_in': 0.1, 'x_caustic_out': 0.15,
                              't_caustic_in': 1000, 't_caustic_out': 2000.,
                              's': 1, 'q': 0.8})
    assert p_1.n_lenses == 1
    assert p_2.n_lenses == 2
    assert p_3.n_lenses == 2


def test_is_finite_source():
    """
    Test if .is_finite_source() works properly for 1L1S
    """
    common = {'t_0': 1.23, 'u_0': 0.123, 't_E': 23.456}
    params_1 = mm.ModelParameters(common)
    params_2 = mm.ModelParameters({'rho': 0.001, **common})
    params_3 = mm.ModelParameters({'t_star': 0.012, **common})

    assert not params_1.is_finite_source()
    assert params_2.is_finite_source()
    assert params_3.is_finite_source()
