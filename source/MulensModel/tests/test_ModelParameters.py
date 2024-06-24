import unittest
import warnings
import pytest
import numpy as np

from astropy import units as u

import MulensModel as mm


class TestModelParameters(unittest.TestCase):
    def test_too_many_parameters_for_init(self):
        """
        Make sure that over-defining parallax fails.
        """
        params = {'t_0': 0, 't_E': 1., 'u_0': 0.1}
        with self.assertRaises(KeyError):
            mm.ModelParameters({**params, 'pi_E': (1., 1.), 'pi_E_N': 1.})
        with self.assertRaises(KeyError):
            mm.ModelParameters({**params, 'pi_E': (1., 1.), 'pi_E_E': 1.})

    def test_wrong_type_of_parameters(self):
        """
        Make sure type of parameters is correct and ranges are also fine
        """
        with self.assertRaises(TypeError):
            mm.ModelParameters({'t_0': 'abc', 'u_0': 1, 't_E': 10.})
        with self.assertRaises(TypeError):
            mm.ModelParameters({'t_0': 123., 'u_0': 1, 't_E': 10.,
                                'shear_G': 0.1, 'alpha': 123.})
        with self.assertRaises(ValueError):
            mm.ModelParameters({'t_0': 123., 'u_0': 1, 't_E': 10., 's': 1.2,
                               'alpha': 34.56, 'q': -0.5})

    def test_init_for_2_sources(self):
        """
        Test different problems for binary source models that should fail.
        """
        with self.assertRaises(Exception):
            mm.ModelParameters({'t_0_1': 1, 'u_0_1': 0.1, 't_E': 10,
                                't_0_2': 10., 'u_0_2': 0.01, 'rho_2': -3.})
        with self.assertRaises(Exception):
            mm.ModelParameters({'t_0_1': 1, 'u_0_1': 0.1, 't_E': 10,
                                'rho_1': -0.1, 't_0_2': 10., 'u_0_2': 0.01})
        with self.assertRaises(KeyError):
            mm.ModelParameters({'t_01': 1, 'u_0_1': 0.1, 't_eff_1': 10,
                                't_0_2': 10., 'u_0_2': 0.01, 't_eff_2': 20.})
        with self.assertRaises(Exception):
            mm.ModelParameters({'t_0_1': 1, 'u_0_2': 0.1, 't_E': 10})


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
    out_2 = "2456141.59300  0.542500    62.6300 "

    assert (out_1 + out_2) == str(params)


def test_repr_no_t_0_par():
    """
    Make sure that t_0_par is printed even if not provided directly.
    """
    t_0 = 2456141.
    u_0 = 0.01
    t_E = 62.63
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E,
                                 'pi_E_E': 0.1, 'pi_E_N': -0.2})

    out_1 = ("    t_0 (HJD)       u_0    t_E (d)    pi_E_N    pi_E_E "
             "t_0_par (HJD) \n")
    out_2 = ("2456141.00000  0.010000    62.6300  -0.20000   0.10000 "
             "2456141.00000 ")

    assert (out_1 + out_2) == str(params)


def test_repr_t_0_par():
    """
    Make sure that t_0_par is printed properly if provided directly.
    """
    t_0 = 2456141.
    u_0 = 0.01
    t_E = 62.63
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E,
                                 'pi_E_E': 0.1, 'pi_E_N': -0.2,
                                 't_0_par': t_0+1})

    out_1 = ("    t_0 (HJD)       u_0    t_E (d)    pi_E_N    pi_E_E "
             "t_0_par (HJD) \n")
    out_2 = ("2456141.00000  0.010000    62.6300  -0.20000   0.10000 "
             "2456142.00000 ")

    assert (out_1 + out_2) == str(params)


def test_repr_t_0_kep():
    """
    Make sure that t_0_kep is printed properly if provided directly.
    """
    t_0 = 2456145.
    u_0 = 0.01
    t_E = 62.63
    s = 1.0
    q = 0.003
    alpha = 30.
    params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,
         'ds_dt': 0.47, 'dalpha_dt': 3.14, 't_0_kep': t_0+1})

    out_1 = ("    t_0 (HJD)       u_0    t_E (d)         s            q "
             "alpha (deg) ds/dt (/yr) dalpha/dt (deg/yr) t_0_kep (HJD) \n")
    out_2 = ("2456145.00000  0.010000    62.6300   1.00000   0.00300000 "
             "   30.00000     0.47000            3.14000 2456146.00000 ")

    assert (out_1 + out_2) == str(params)


def test_positive_t_E():
    """
    Check if t_E is positive when t_eff is given, even if u_0 is negative.
    """
    t_0 = 10205.1
    u_0 = -0.50
    t_eff = 12.5
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_eff': t_eff})

    assert params.t_E >= 0.
    assert params.t_E == params.t_eff / abs(params.u_0)


def test_q_gt_1_is_good():
    """
    Check if the magnification is reproduced by transforming q -> 1/q and
    alpha -> alpha +/- 180, including for q > 1. See issue #84.
    """
    t_0 = 3583.
    u_0 = 0.3
    t_E = 12.
    s = 1.65
    q = 0.25
    alpha = 339.0
    rho = 0.001

    planet = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q,
                       'alpha': alpha, 'rho': rho})
    planet_2 = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': 1/q,
                         'alpha': alpha+180., 'rho': rho})
    planet_3 = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': 1/q,
                         'alpha': alpha-180., 'rho': rho})
    list_of_methods = [3588., 'VBBL', 3594., 'hexadecapole', 3598.0]
    planet.set_magnification_methods(list_of_methods)
    planet_2.set_magnification_methods(list_of_methods)
    planet_3.set_magnification_methods(list_of_methods)
    t_checks = [3580, 3589, 3590, 3592, 3593, 3595]
    magnifications_1 = planet.get_magnification(time=t_checks)
    magnifications_2 = planet_2.get_magnification(time=t_checks)
    magnifications_3 = planet_3.get_magnification(time=t_checks)

    assert max(magnifications_1 - magnifications_2) < 1e-10
    assert max(magnifications_1 - magnifications_3) < 1e-10


def test_q_gt_1_is_smooth():
    """
    Check that there is a smooth transition between q = 0.97, 0.99 and 1.01.
    In this case, a trajectory with caustic approaching is adopted.
    """
    t_0 = 3583.
    u_0 = 0.3
    t_E = 12.
    s = 2.18
    q = 0.99
    alpha = 310.
    rho = 0.001
    planet = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q,
                       'alpha': alpha, 'rho': rho})
    q_min = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q-0.02,
                      'alpha': alpha, 'rho': rho})
    q_max = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q+0.02,
                      'alpha': alpha, 'rho': rho})

    planet.set_magnification_methods([3580., 'VBBL', 3595.])
    q_min.set_magnification_methods([3580., 'VBBL', 3595.])
    q_max.set_magnification_methods([3580., 'VBBL', 3595.])
    t_checks = [3571, 3583, 3585.5, 3586, 3586.5, 3592.5]
    magnification = planet.get_magnification(time=t_checks)
    diff_min = magnification - q_min.get_magnification(time=t_checks)
    diff_max = magnification - q_max.get_magnification(time=t_checks)
    limits = np.array([0.01, 0.01, 0.01, 0.56, 0.01, 0.03])

    assert np.all(abs(diff_min) < limits)
    assert np.all(abs(diff_max) < limits)


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
        t_0 = 2450000.
        u_0 = 0.1
        t_E = 20. * u.day
        rho = 0.001
        t_star = t_E * rho
        with self.assertRaises(KeyError):
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


def test_single_lens_convergence_K_shear_G():
    """
    Test single lens with convergence_K and shear_G in intialized
    """
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    convergence_K = 0.1
    shear_G = complex(0.1, 0.2)
    alpha = 20.
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'convergence_K': convergence_K,
        'shear_G': shear_G, 'alpha': alpha})

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E.value)
    np.testing.assert_almost_equal(params.convergence_K, convergence_K)
    np.testing.assert_almost_equal(params.shear_G.real, shear_G.real)

    convergence_K *= 2.
    shear_G *= 2.
    params.convergence_K = convergence_K
    params.shear_G = shear_G
    np.testing.assert_almost_equal(params.convergence_K, convergence_K)
    np.testing.assert_almost_equal(params.shear_G.real, shear_G.real)


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


def test_single_lens_with_mass_sheet():
    """
    Test if Chang-Refsdal microlensing parameters are properly defined.
    """
    basic = {'t_0': 1000., 'u_0': 0.1, 't_E': 20.}
    G = complex(-0.1, -0.2)
    K = -0.1

    _ = mm.ModelParameters({**basic})
    _ = mm.ModelParameters({**basic, 'shear_G': G, 'alpha': 123.})
    _ = mm.ModelParameters({**basic, 'convergence_K': K})
    _ = mm.ModelParameters(
        {**basic, 'shear_G': G, 'convergence_K': K, 'alpha': 123.})


class TestParameters(unittest.TestCase):
    def test_failing_single_lens_with_mass_sheet(self):
        """
        Test if Chang-Refsdal microlensing fails when it's expected to fail.
        """
        basic = {'t_0': 1000., 'u_0': 0.1, 't_E': 20.}
        G = complex(-0.1, -0.2)
        K = -0.1
        alpha = 123.

        # Cases with missing one of PSPL parameters:
        dict_K = {'convergence_K': K}
        dict_G = {'shear_G': G}
        for mass_sheet in [dict_G, dict_K, {**dict_G, **dict_K}]:
            with self.assertRaises(KeyError):
                mm.ModelParameters({'t_0': 1000., 'u_0': 0.1, **mass_sheet})
            with self.assertRaises(KeyError):
                mm.ModelParameters({'t_E': 20., 'u_0': 0.1, **mass_sheet})
            with self.assertRaises(KeyError):
                mm.ModelParameters({'t_0': 1000., 't_E': 20., **mass_sheet})

        # Cases below have too many parameters:
        with self.assertRaises(KeyError):
            mm.ModelParameters({**basic, 'alpha': alpha})
        with self.assertRaises(KeyError):
            mm.ModelParameters({**basic, 'dalpha_dt': 0.01})
        with self.assertRaises(KeyError):
            mm.ModelParameters(
                {**basic, 'convergence_K': K, 'alpha': alpha})
        with self.assertRaises(KeyError):
            mm.ModelParameters({**basic, 'convergence_K': K,
                                'alpha': alpha, 'dalpha_dt': -0.3})
        with self.assertRaises(KeyError):
            mm.ModelParameters(
                {**basic, 'shear_G': G, 'convergence_K': K, 'alpha': alpha,
                 'dalpha_dt': -0.3})

        # The case below is missing alpha:
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters({**basic, 'shear_G': G})

        # No access to mass sheet without it's parameters
        parameters = mm.ModelParameters(basic)
        with self.assertRaises(KeyError):
            parameters.convergence_K = 0.
        with self.assertRaises(KeyError):
            parameters.convergence_K
        with self.assertRaises(KeyError):
            parameters.shear_G = 0.


xallarap_parameters = {
    't_0': 2., 't_E': 9., 'u_0': 0.1, 'xi_period': 12.345,
    'xi_semimajor_axis': 0.54321, 'xi_Omega_node': 0.123,
    'xi_inclination': 9.8765, 'xi_argument_of_latitude_reference': 24.68,
    'xi_eccentricity': 0.5, 'xi_omega_periapsis': 12.3456, 't_0_xi': 1.}


def setup_xallarap(key):
    """
    Setup for xallarap tests.
    """
    model = mm.ModelParameters({**xallarap_parameters})
    return (model, xallarap_parameters[key])


tested_keys_1 = ['xi_period', 'xi_semimajor_axis', 'xi_Omega_node',
                 'xi_inclination', 'xi_argument_of_latitude_reference',
                 'xi_eccentricity', 'xi_omega_periapsis']
tested_keys_2 = tested_keys_1 + ['t_0_xi']


@pytest.mark.parametrize("key", tested_keys_2)
def test_xallarap_set_value_1(key):
    """
    Check if xallarap settings of xallarap are correctly changed via dictionary
    """
    (model, value) = setup_xallarap(key)
    new_value = 2.1234 * value
    model.parameters[key] = new_value
    assert model.parameters[key] == new_value
    assert getattr(model, key) == new_value


@pytest.mark.parametrize("key", tested_keys_2)
def test_xallarap_set_value_2(key):
    """
    Check if xallarap settings are correctly changed via attribute
    """
    (model, value) = setup_xallarap(key)
    new_value = 1.3579 * value
    setattr(model, key, new_value)
    assert model.parameters[key] == new_value
    assert getattr(model, key) == new_value


class TestXallarapErrors(unittest.TestCase):
    def test_missing_xallarap_parameters(self):
        """Make sure that the required xallarap parameters are all defined"""
        for parameter in tested_keys_1:
            parameters = {**xallarap_parameters}
            parameters.pop(parameter)
            with self.assertRaises(KeyError):
                mm.ModelParameters(parameters)

    def test_negative_period(self):
        """
        Make sure that period is positive
        """
        key = 'xi_period'
        (model, value) = setup_xallarap(key)
        new_value = -3.14 * value
        with self.assertRaises(ValueError):
            setattr(model, key, new_value)

    def test_negative_semimajor_axis(self):
        """
        Make sure that period is positive
        """
        key = 'xi_semimajor_axis'
        (model, value) = setup_xallarap(key)
        new_value = -3.14 * value
        with self.assertRaises(ValueError):
            setattr(model, key, new_value)

    def test_xallarap_and_binary_source(self):
        """
        Confirm that binary source and xallrap cannot be defined in
        the same model
        """
        parameters = {
            't_0_1': 0, 'u_0_1': 1, 't_0_2': 5, 'u_0_2': 0.1, 't_E': 9,
            'xi_period': 12., 'xi_semimajor_axis': 0.5, 'xi_Omega_node': 10.,
            'xi_inclination': 50., 'xi_argument_of_latitude_reference': 200.}
        with self.assertRaises(NotImplementedError):
            mm.ModelParameters(parameters)

    def test_xallarap_and_Cassan08(self):
        """
        Confirm that xallrap and binary lens in Cassan 2008 parameterization
        cannot be defined jointly
        """
        parameters = {
            's': 1, 'q': 0.8, 'x_caustic_in': 0.1, 'x_caustic_out': 0.15,
            't_caustic_in': 1000, 't_caustic_out': 2000.,
            'xi_period': 12., 'xi_semimajor_axis': 0.5, 'xi_Omega_node': 10.,
            'xi_inclination': 50., 'xi_argument_of_latitude_reference': 200.}
        with self.assertRaises(NotImplementedError):
            mm.ModelParameters(parameters)

    def test_negative_source_mass_ratio_1(self):
        """
        q_source must be positive in __init__()
        """
        parameters = {**xallarap_parameters, 'q_source': -0.12345}
        with self.assertRaises(ValueError):
            _ = mm.ModelParameters(parameters)

    def test_negative_source_mass_ratio_2(self):
        """
        q_source must be positive
        """
        parameters = {**xallarap_parameters, 'q_source': 0.12345}
        model = mm.ModelParameters(parameters)
        with self.assertRaises(ValueError):
            setattr(model, 'q_source', -0.12345)

    def test_overdefined_source_size(self):
        """
        overdefine first sourece size
        """
        parameters = {**xallarap_parameters,
                      'rho_1': 0.1, 't_star_1': 0.1, 'q_source': 1.0}
        with self.assertRaises(ValueError):
            mm.ModelParameters(parameters)

    def test_mixed_binary_source_and_xallarap(self):
        """
        Xallarap cannot be combined with just
        one stardard binary source parameter.
        """
        parameters = {**xallarap_parameters, 'u_0_2': 0.123}
        with self.assertRaises(ValueError):
            mm.ModelParameters(parameters)

    def test_no_q_source_but_with_rho_2(self):
        """
        xallarap model without q_source cannot have rho_2
        """
        parameters = {**xallarap_parameters, 'rho_2': 0.1}
        with self.assertRaises(KeyError):
            mm.ModelParameters(parameters)

    def test_PSPL_and_q_source(self):
        """
        Make sure one cannot provide only PSPL parameters and q_source.
        """
        with self.assertRaises(KeyError):
            mm.ModelParameters({'t_0': 1, 'u_0': 2, 't_E': 3, 'q_source': 1})


def test_print_xallarap():
    """
    Test if printing of printing of xallarap model works as expected.
    """
    model = mm.ModelParameters(xallarap_parameters)
    expected = (
        "    t_0 (HJD)       u_0    t_E (d) xallarap period (d) xallarap "
        "semimajor axis xallarap inclination (deg) xallarap Omega node (deg) "
        "xallarap argument of latitude reference (deg) xallarap eccentricity "
        "xallarap omega periapsis (deg)  t_0_xi (HJD) "
        "\n      2.00000  0.100000     9.0000             12.3450            "
        "    0.543210                    9.87650                   0.12300   "
        "                                   24.68000              0.500000   "
        "                    12.34560       1.00000 "
        "\nxallarap reference position: (0.2673, 0.0582)"
        )
    assert model.__repr__() == expected


def test_print_xallarap_with_q_source():
    """
    Test if printing of printing of xallarap model with q_source works
    as expected. Most stuff was tested in test_print_xallarap(), so we
    check only the parts that are important here.
    """
    parameters = {**xallarap_parameters, 'q_source': 0.12345}
    model = mm.ModelParameters(parameters)
    lines = model.__repr__().split("\n")
    assert lines[0][-27:] == "    q_source  t_0_xi (HJD) "
    assert lines[1][-27:] == "  0.12345000       1.00000 "
    assert lines[2] == "xallarap reference position: (0.2673, 0.0582)"


@pytest.mark.parametrize(
    "parameter",
    ['xi_Omega_node', 'xi_inclination', 'xi_argument_of_latitude_reference'])
@pytest.mark.parametrize("value", [-361., 541.])
def test_warnings_for_xallarap_angles(parameter, value):
    """
    Check if xallarap angles in somehow strange range give warning
    """
    parameters = {
        't_0': 0, 't_E': 9., 'u_0': 0.1, 'xi_period': 12.345,
        'xi_semimajor_axis': 0.54321, 'xi_Omega_node': 100.,
        'xi_inclination': 50., 'xi_argument_of_latitude_reference': 200.,
        't_0_xi': 1.}
    parameters[parameter] = value

    with warnings.catch_warnings(record=True) as warnings_:
        warnings.simplefilter("always")
        mm.ModelParameters(parameters)
        assert len(warnings_) == 1
        assert issubclass(warnings_[0].category, RuntimeWarning)


def test_is_xallarap_1():
    """
    Make sure that .is_xallarap() returns True, when it should.
    """
    parameters = {
        't_0': 0, 't_E': 9., 'u_0': 0.1, 'xi_period': 12.345,
        'xi_semimajor_axis': 0.54321, 'xi_Omega_node': 100.,
        'xi_inclination': 50., 'xi_argument_of_latitude_reference': 200.,
        't_0_xi': 1.}
    model_params = mm.ModelParameters(parameters)
    assert model_params.is_xallarap


def test_is_xallarap_2():
    """
    Check that is_xallarap() returns False for a (static) binary source model.
    """
    parameters = {'t_0_1': 0, 'u_0_1': 1, 't_0_2': 5, 'u_0_2': 0.1, 't_E': 9}
    model_params = mm.ModelParameters(parameters)
    assert not model_params.is_xallarap


def test_xallarap_n_sources():
    """
    Make sure that number of sources in xallarap models is properly calculated
    """
    parameters = {**xallarap_parameters}
    model_1S = mm.ModelParameters(parameters)
    assert model_1S.n_sources == 1

    parameters['q_source'] = 1.
    model_2S = mm.ModelParameters(parameters)
    assert model_2S.n_sources == 2

    parameters['rho_1'] = 0.1
    parameters['rho_2'] = 0.2
    model_3 = mm.ModelParameters(parameters)
    assert model_3.n_sources == 2

    parameters = {**xallarap_parameters, 't_star_1': 2.}
    model_4 = mm.ModelParameters(parameters)
    assert model_4.n_sources == 1


def _test_2S1L_xallarap_individual_source_parameters(xi_u):
    """
    Make sure that parameters of both sources are properly set.
    Most importantly, xi_u is shifted by 180 deg and xi_a is scaled by q_source.
    """
    q_source = 1.23456
    parameters_1st = {**xallarap_parameters}
    parameters_1st['xi_argument_of_latitude_reference'] = xi_u

    parameters_2nd = {**parameters_1st}
    parameters_2nd['xi_semimajor_axis'] /= q_source
    if xi_u < 180:
        parameters_2nd['xi_argument_of_latitude_reference'] += 180.
    else:
        parameters_2nd['xi_argument_of_latitude_reference'] -= 180.

    parameters = {'q_source': q_source, **parameters_1st}
    model = mm.ModelParameters(parameters)
    check_1st = model.source_1_parameters.as_dict()
    check_1st['t_E'] = check_1st['t_E'].value
    check_2nd = model.source_2_parameters.as_dict()
    check_2nd['t_E'] = check_2nd['t_E'].value

    assert check_1st == parameters_1st
    assert check_2nd == parameters_2nd


def test_2S1L_xallarap_individual_source_parameters_1():
    """
    Make sure xi_u is increased by 180 for small input.
    """
    _test_2S1L_xallarap_individual_source_parameters(xi_u=8.642)


def test_2S1L_xallarap_individual_source_parameters_2():
    """
    Make sure xi_u is increased by 180 for large input.
    """
    _test_2S1L_xallarap_individual_source_parameters(xi_u=234.567)


tested_keys_3 = tested_keys_2 + ['q_source']


@pytest.mark.parametrize("key", tested_keys_3)
def test_changes_of_xallrap_parameters_for_both_sources(key):
    """
    Make sure that chainging a xallarap parameter in a binary source event
    with binary sources model properly changes parameters of each parameter.
    For q_source make sure that it actually it's not passed to the parameters
    of each source.
    """
    q_source = 1.23456
    factor = 1.1
    parameters = {'q_source': q_source, **xallarap_parameters}
    model = mm.ModelParameters(parameters)
    old_value = getattr(model, key)
    new_value = factor * old_value
    setattr(model, key, new_value)

    assert getattr(model, key) == new_value

    if key == 'q_source':
        assert 'q_source' not in model.source_1_parameters.parameters
        assert 'q_source' not in model.source_2_parameters.parameters
        key_a = 'xi_semimajor_axis'
        xi_a = xallarap_parameters[key_a]
        assert model.source_1_parameters.parameters[key_a] == xi_a
        assert model.source_2_parameters.parameters[key_a] == xi_a / new_value
        return

    assert getattr(model.source_1_parameters, key) == new_value

    new_value_2 = new_value
    if key == 'xi_argument_of_latitude_reference':
        new_value_2 += 180.
    elif key == 'xi_semimajor_axis':
        new_value_2 /= q_source

    assert getattr(model.source_2_parameters, key) == new_value_2


def test_reference_position():
    """
    Make sure that for xallarap model bith sources have the same
    reference position.
    """
    parameters = {**xallarap_parameters, 'q_source': 0.12345}
    model = mm.ModelParameters(parameters)
    text_1 = model.source_1_parameters.__repr__().split("\n")[-1]
    text_2 = model.source_2_parameters.__repr__().split("\n")[-1]
    assert text_1 == text_2
