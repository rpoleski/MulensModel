import unittest
import pytest
import numpy as np

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
                                'shear_G': 0.1, 'alpha': 303.})
        with self.assertRaises(ValueError):
            mm.ModelParameters({'t_0': 123., 'u_0': 1, 't_E': 10., 's': 1.2,
                               'alpha': 214.56, 'q': -0.5})

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
    t_E = 62.63
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E)


def test_repr_parameters():
    """
    Test str(ModelParameters) or __repr__(ModelParameters)
    """
    t_0 = 2456141.593
    u_0 = 0.5425
    t_E = 62.63
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


class TestT0X(unittest.TestCase):
    """
    Test that t_0_par and t_0_kep are correctly set even if not provided.
    """
    def setUp(self):
        t_0 = 2456145.
        u_0 = 0.01
        t_E = 62.63
        s = 1.0
        q = 0.003
        alpha = 30.
        self.params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q,
                       'alpha': alpha, 'ds_dt': 0.47, 'dalpha_dt': 3.14,
                       'pi_E_E': 0.1, 'pi_E_N': -0.2}

    def test_basic(self):
        model_params = mm.ModelParameters(self.params)
        np.testing.assert_almost_equal(model_params.t_0_par, self.params['t_0'])
        np.testing.assert_almost_equal(model_params.t_0_kep, self.params['t_0'])

    def test_N_sources(self):
        params = self.params.copy()
        params['t_0_1'] = self.params['t_0'] - 1.
        params['u_0_1'] = self.params['u_0'] * 2
        params['t_0_2'] = self.params['t_0'] + 1.
        params['u_0_2'] = 0.1
        params.pop('t_0')
        params.pop('u_0')

        assert 't_0' not in params.keys()
        assert 'u_0' not in params.keys()

        model_params = mm.ModelParameters(params)

        np.testing.assert_almost_equal(model_params.t_0_par, self.params['t_0'] - 1.)
        np.testing.assert_almost_equal(model_params.t_0_kep, self.params['t_0'] - 1.)

    def test_t_0_par_set(self):
        params = self.params.copy()
        params['t_0_par'] = self.params['t_0'] + 3.

        assert 't_0_par' in params.keys()
        assert 't_0_kep' not in params.keys()

        model_params = mm.ModelParameters(params)
        np.testing.assert_almost_equal(model_params.t_0_par, self.params['t_0'] + 3.)
        np.testing.assert_almost_equal(model_params.t_0_kep, self.params['t_0'] + 3.)

    def test_t_0_kep_set(self):
        params = self.params.copy()
        params['t_0_kep'] = self.params['t_0'] + 5.

        assert 't_0_par' not in params.keys()
        assert 't_0_kep' in params.keys()

        model_params = mm.ModelParameters(params)
        np.testing.assert_almost_equal(model_params.t_0_par, self.params['t_0'] + 5.)
        np.testing.assert_almost_equal(model_params.t_0_kep, self.params['t_0'] + 5.)


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
    alpha = 159.0
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

    assert max(magnifications_1 - magnifications_2) < 1e-7
    assert max(magnifications_2 - magnifications_3) < 1e-10


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
    alpha = 130.
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
    t_E_2 = t_E_1
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
        t_E = 20.
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
    t_E = 62.63
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    params.as_dict().update({'rho': 0.001})

    assert len(params.parameters.keys()) == 4


def test_orbital_motion_1():
    """basic tests of orbital motion"""
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 210.}
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
    assert static.get_alpha(epoch_1) == dict_static['alpha']

    # Test get_s() and get_alpha() for orbital motion case.
    np.testing.assert_almost_equal(motion.get_alpha(epoch_1), 209.5)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_2), 210.)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_3), 210.5)
    np.testing.assert_almost_equal(motion.get_s(epoch_1), 1.2295)
    np.testing.assert_almost_equal(motion.get_s(epoch_2), 1.2345)
    np.testing.assert_almost_equal(motion.get_s(epoch_3), 1.2395)

    # Test arguments as list or array.
    np.testing.assert_almost_equal(
        motion.get_alpha([epoch_1, epoch_2, epoch_3]),
        [209.5, 210., 210.5])
    np.testing.assert_almost_equal(
        motion.get_alpha(np.array([epoch_1, epoch_2, epoch_3])),
        [209.5, 210., 210.5])
    np.testing.assert_almost_equal(
        motion.get_s([epoch_1, epoch_2, epoch_3]),
        [1.2295, 1.2345, 1.2395])
    np.testing.assert_almost_equal(
        motion.get_s(np.array([epoch_1, epoch_2, epoch_3])),
        [1.2295, 1.2345, 1.2395])

    assert motion.alpha == 210.
    assert motion.s == 1.2345


def test_t_0_kep():
    """test reference epoch for orbital motion"""
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 210.}
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
    np.testing.assert_almost_equal(motion.get_alpha(epoch_1), 209.5)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_2), 210.)
    np.testing.assert_almost_equal(motion.get_s(epoch_1), 1.2295)
    np.testing.assert_almost_equal(motion.get_s(epoch_2), 1.2345)

    # Test motion_2.
    np.testing.assert_almost_equal(motion_2.get_alpha(epoch_1), 210.)
    np.testing.assert_almost_equal(motion_2.get_alpha(epoch_2), 210.5)
    np.testing.assert_almost_equal(motion_2.get_s(epoch_1), 1.2345)
    np.testing.assert_almost_equal(motion_2.get_s(epoch_2), 1.2395)


def setup_keplerian(dict_to_add):
    """
    Setup dictionary for tests of circular keplerian motion
    """
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 30., 'rho': 0.001}

    return {**dict_static, **dict_to_add}


class test_Keplerian(unittest.TestCase):
    def test_keplerian_both_inputs(self):
        """s_z and ds_z_dt are given with ds_dt and dalpha_dt"""
        dict_2 = setup_keplerian({'ds_dt': 0.1, 'dalpha_dt': 10.,
                                  's_z': 0.1, 'ds_z_dt': 1.9})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_2)

    def test_keplerian_no_ds_dt(self):
        """fails if s_z and ds_z_dt are given with dalpha_dt"""
        dict_3 = setup_keplerian({'dalpha_dt': 10., 's_z': 0.1,
                                  'ds_z_dt': 1.9})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_3)

    def test_keplerian_only_z(self):
        """fails if s_z and ds_z_dt are given only"""
        dict_4 = setup_keplerian({'s_z': 0.1, 'ds_z_dt': 1.9})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_4)


def test_is_not_keplerian():
    """
    Make sure that is_keplerian() works properly for static model.
    """
    params_static = setup_keplerian({})
    model_static = mm.ModelParameters(params_static)
    assert model_static.is_static()
    assert not model_static.is_keplerian()


def test_keplerian_motion_s_z_only():
    """tests if only s_z is given"""
    dict_1 = setup_keplerian({'ds_dt': 0.1, 'dalpha_dt': 10., 's_z': 0.1})
    keplerian = mm.ModelParameters(dict_1)
    assert not keplerian.is_static()
    assert keplerian.is_keplerian()


def test_keplerian_motion_ds_z_dt_only():
    """tests if only ds_z_dt is given"""
    dict_1 = setup_keplerian({'ds_dt': 0.1, 'dalpha_dt': 10., 'ds_z_dt': 1.9})
    keplerian = mm.ModelParameters(dict_1)
    assert not keplerian.is_static()
    assert keplerian.is_keplerian()


def test_circular_motion_setting_s_z():
    """
    Check if s_z is properly set.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 's_z': 0.1})
    model = mm.ModelParameters(params)
    assert model.s_z == 0.1


def test_circular_motion_setting_ds_z_dt():
    """
    Check if ds_z_dt is properly set.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 'ds_z_dt': 0.1})
    model = mm.ModelParameters(params)
    assert model.ds_z_dt == 0.1


def test_calculation_ds_z_dt_for_circular_motion_1():
    """
    Set s_z for circular motion and test if ds_z_dt is calculated properly.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 's_z': 0.1})
    model = mm.ModelParameters(params)
    np.testing.assert_almost_equal(model.ds_z_dt, -2.469)


def test_calculation_s_z_for_circular_motion_1():
    """
    Set ds_z_dt for circular motion and test if s_z is calculated properly.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 'ds_z_dt': 0.1})
    model = mm.ModelParameters(params)
    np.testing.assert_almost_equal(model.s_z, -2.469)


def test_calculation_ds_z_dt_for_circular_motion_2():
    """
    Set s_z for circular motion and test if ds_z_dt is calculated properly.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 's_z': 0.01})
    model = mm.ModelParameters(params)
    model.s_z = 0.1
    np.testing.assert_almost_equal(model.ds_z_dt, -2.469)


def test_calculation_s_z_for_circular_motion_2():
    """
    Set ds_z_dt for circular motion and test if s_z is calculated properly.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 'ds_z_dt': 1.0})
    model = mm.ModelParameters(params)
    model.ds_z_dt = 0.1
    np.testing.assert_almost_equal(model.s_z, -2.469)


def test_calculation_s_z_for_circular_motion_3():
    """
    Set ds_z_dt for circular motion and test if s_z is calculated properly.
    This time access gamma_parallel in the meantime.
    This test was failing at some point.
    """
    dict_params = setup_orbital_motion_gammas(
        {'dalpha_dt': 10, 'ds_dt': 0., 'ds_z_dt': 10*2**-0.5*(np.pi/180.)})
    model = mm.ModelParameters(dict_params)
    model.gamma_parallel
    np.testing.assert_almost_equal(model.s_z, 0.)


def test_access_to_ds_z_dt():
    """
    Make sure that accessing ds_z_dt doesn't change s_z.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 's_z': 0.01})
    model = mm.ModelParameters(params)
    model.s_z = 0.1
    model.ds_z_dt
    assert model.s_z == 0.1


def test_access_to_s_z():
    """
    Make sure that accessing s_z doesn't change ds_z_dt.
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 'ds_z_dt': 0.01})
    model = mm.ModelParameters(params)
    model.ds_z_dt = 0.1
    model.s_z
    assert model.ds_z_dt == 0.1


def test_access_to_s_z_2():
    """
    Make sure that accessing s_z doesn't fix it
    (i.e., it changes when ds_z_dt changes).
    """
    params = setup_keplerian({'ds_dt': 0.2, 'dalpha_dt': 10., 'ds_z_dt': 0.01})
    model = mm.ModelParameters(params)
    np.testing.assert_almost_equal(model.s_z, -24.69)
    model.ds_z_dt = 0.1
    np.testing.assert_almost_equal(model.s_z, -2.469)


def setup_orbital_motion_gammas(dict_to_add):
    """
    Setup dictionary for tests of gammas in orbital motion
    """
    dict_pars = {'t_0': 2457123.456, 'u_0': 0.0345, 't_E': 30.00, 's': 1.5,
                 'q': 0.987, 'alpha': 12.345}

    return {**dict_pars, **dict_to_add}


def test_orbital_motion_gammas():
    """test .gamma_parallel .gamma_perp .gamma"""

    dict_params = setup_orbital_motion_gammas({'ds_dt': 0.5, 'dalpha_dt': 50.})
    params = mm.ModelParameters(dict_params)

    np.testing.assert_almost_equal(params.gamma_parallel, 0.333333333)
    np.testing.assert_almost_equal(params.gamma_perp, -0.872664626)
    np.testing.assert_almost_equal(params.gamma, 0.934159869)


def test_orbital_motion_gamma_z():
    """test .gamma_z .gamma"""
    dict_pars = setup_orbital_motion_gammas({'ds_dt': 0.5, 'dalpha_dt': 50.,
                                             's_z': 0.5})
    params = mm.ModelParameters(dict_pars)

    np.testing.assert_almost_equal(params.gamma_z, -1.)
    np.testing.assert_almost_equal(params.gamma, 1.368449729)


def test_lens_orbital_parameters_circular_1():
    """
    Check if parameters of face-on circular orbital motion
    are properly calculated.
    """
    dict_params = setup_orbital_motion_gammas(
        {'dalpha_dt': -36./2**.5, 'ds_dt': 0., 'ds_z_dt': 36./2**.5*(np.pi/180.)*1.5})
    dict_params['alpha'] = 90.

    parameters = mm.ModelParameters(dict_params)

    np.testing.assert_almost_equal(parameters.lens_semimajor_axis, 1.5)
    np.testing.assert_almost_equal(parameters.lens_period, 10.)
    np.testing.assert_almost_equal(parameters.lens_inclination, 45.)
    np.testing.assert_almost_equal(parameters.lens_Omega_node, 0.)
    np.testing.assert_almost_equal(parameters.lens_argument_of_latitude_reference, 0.)


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
                              's': 10, 'q': 0.5, 'alpha': 280.})
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
    t_E = 62.63
    convergence_K = 0.1
    shear_G = complex(0.1, 0.2)
    alpha = 200.
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'convergence_K': convergence_K,
        'shear_G': shear_G, 'alpha': alpha})

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E)
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
    _ = mm.ModelParameters({**basic, 'shear_G': G, 'alpha': 303.})
    _ = mm.ModelParameters({**basic, 'convergence_K': K})
    _ = mm.ModelParameters(
        {**basic, 'shear_G': G, 'convergence_K': K, 'alpha': 303.})


class TestParameters(unittest.TestCase):
    def test_failing_single_lens_with_mass_sheet(self):
        """
        Test if Chang-Refsdal microlensing fails when it's expected to fail.
        """
        basic = {'t_0': 1000., 'u_0': 0.1, 't_E': 20.}
        G = complex(-0.1, -0.2)
        K = -0.1
        alpha = 303.

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
        with self.assertRaises(AttributeError):
            parameters.convergence_K = 0.
        with self.assertRaises(AttributeError):
            parameters.convergence_K
        with self.assertRaises(AttributeError):
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
                print("Test failed (i.e. KeyError was not raised) for ",
                      parameter)

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
        with self.assertRaises(KeyError):
            mm.ModelParameters(parameters)

    def test_rho_with_source_unspecified(self):
        parameters = {**xallarap_parameters, 'q_source': 1.3, 'rho': 0.1}
        with self.assertRaises(KeyError):
            mm.ModelParameters(parameters)

    def test_mixed_binary_source_and_xallarap(self):
        """
        Xallarap cannot be combined with just
        one stardard binary source parameter.
        """
        parameters = {**xallarap_parameters, 'u_0_2': 0.123}
        with self.assertRaises(KeyError):
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
    Most importantly, xi_u is shifted by 180deg and xi_a is scaled by q_source.
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
    check_2nd = model.source_2_parameters.as_dict()

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
    Make sure that for xallarap model, both sources have the same
    reference position.
    """
    parameters = {**xallarap_parameters, 'q_source': 0.12345}
    model = mm.ModelParameters(parameters)
    text_1 = model.source_1_parameters.__repr__().split("\n")[-1]
    text_2 = model.source_2_parameters.__repr__().split("\n")[-1]
    assert text_1 == text_2


def test_change_of_xallarap_reference_position_1():
    """
    Make sure that changing xallarap parameters changes the reference position.
    """
    parameters = setup_xallarap('t_0')[0]
    reference_1 = parameters.xallarap_reference_position
    parameters.xi_omega_periapsis = 45.
    reference_2 = parameters.xallarap_reference_position
    assert np.all(reference_1 != reference_2)


def test_change_of_xallarap_reference_position_2():
    """
    Make sure that the xallarap reference position is based on the current parameters.
    """
    new_value = 45.
    parameters_1 = mm.ModelParameters({**xallarap_parameters})
    parameters_1.xi_omega_periapsis = new_value
    reference_1 = parameters_1.xallarap_reference_position

    changed_values = {**xallarap_parameters, 'xi_omega_periapsis': new_value}
    parameters_2 = mm.ModelParameters({**changed_values})
    reference_2 = parameters_2.xallarap_reference_position

    assert np.all(reference_1 == reference_2)


class Test1L3SModels(unittest.TestCase):

    def setUp(self):
        self.t_0 = [0., 5., 2.]
        self.u_0 = [1., 0.1, 0.3]
        self.t_E = 9.
        rho_2 = 0.001
        t_star_3 = 0.02
        self.rho = [0., rho_2, t_star_3 / self.t_E]
        self.t_star = np.array(self.rho) * self.t_E

        parameters = {'t_0_1': self.t_0[0], 'u_0_1': self.u_0[0],
                      't_0_2': self.t_0[1], 'u_0_2': self.u_0[1],
                      'rho_2': rho_2,
                      't_0_3': self.t_0[2], 'u_0_3': self.u_0[2],
                      't_star_3': t_star_3,
                      't_E': self.t_E}
        self.model_params = mm.ModelParameters(parameters)

    # Basic 1L3S Tests
    def test_basic_1L3S_n_sources(self):
        """
        Test that we can define Triple source models.
        """
        assert self.model_params.n_sources == 3

    def test_basic_1L3S_t_0(self):
        """
        Test that we can access attributes of Triple source models.
        """
        assert (self.model_params.t_0_3 == self.t_0[2])
        for i in range(3):
            assert (self.model_params.__getattr__(
                't_0_{0}'.format(i+1)) == self.t_0[i])

        assert (self.model_params.source_3_parameters.t_0 == self.t_0[2])
        for i in range(3):
            assert (self.model_params.__getattr__(
                'source_{0}_parameters'.format(i+1)).t_0 == self.t_0[i])

        np.testing.assert_almost_equal(self.model_params.t_0_1, self.t_0[0])
        np.testing.assert_almost_equal(self.model_params.t_0_2, self.t_0[1])
        np.testing.assert_almost_equal(self.model_params.t_0_3, self.t_0[2])

    def test_basic_1L3S_u_0(self):
        """
        Test that we can access attributes of Triple source models.
        """
        for i in range(3):
            assert (self.model_params.__getattr__(
                'u_0_{0}'.format(i+1)) == self.u_0[i])

        for i in range(3):
            assert (self.model_params.__getattr__(
                'source_{0}_parameters'.format(i+1)).u_0 == self.u_0[i])

        np.testing.assert_almost_equal(self.model_params.u_0_1, self.u_0[0])
        np.testing.assert_almost_equal(self.model_params.u_0_2, self.u_0[1])
        np.testing.assert_almost_equal(self.model_params.u_0_3, self.u_0[2])

    def test_basic_1L3S_rho(self):
        """
        Test that we can access attributes of Triple source models.
        """
        for i in range(1, 3):
            assert (self.model_params.__getattr__(
                'rho_{0}'.format(i+1)) == self.rho[i])

        for i in range(1, 3):
            assert (self.model_params.__getattr__(
                'source_{0}_parameters'.format(i+1)).rho == self.rho[i])

        with np.testing.assert_raises(AttributeError):
            self.model_params.rho_1

        np.testing.assert_almost_equal(self.model_params.rho_2, self.rho[1])
        np.testing.assert_almost_equal(self.model_params.rho_3, self.rho[2])

    def test_basic_1L3S_t_star(self):
        """
        Test that we can access attributes of Triple source models.
        """
        for i in range(1, 3):
            assert (self.model_params.__getattr__(
                't_star_{0}'.format(i+1)) == self.t_star[i])

        for i in range(1, 3):
            assert (self.model_params.__getattr__(
                'source_{0}_parameters'.format(i+1)).t_star == self.t_star[i])

        with np.testing.assert_raises(AttributeError):
            _ = self.model_params.t_star_1

        np.testing.assert_almost_equal(self.model_params.t_star_2, self.t_star[1])
        np.testing.assert_almost_equal(self.model_params.t_star_3, self.t_star[2])

    def test_teff2u0(self):
        teff = np.array(self.u_0) * self.t_E
        params = {}
        for i in range(3):
            params['t_0_{0}'.format(i+1)] = self.t_0[i]
            params['t_eff_{0}'.format(i+1)] = teff[i]

        params['t_E'] = self.t_E
        model = mm.ModelParameters(params)
        for i in range(3):
            u_0 = model.__getattr__('u_0_{0}'.format(i+1))
            np.testing.assert_almost_equal(u_0, self.u_0[i])

        np.testing.assert_almost_equal(model.u_0_1, self.u_0[0])
        np.testing.assert_almost_equal(model.u_0_2, self.u_0[1])
        np.testing.assert_almost_equal(model.u_0_3, self.u_0[2])


class Test1LNSModels(unittest.TestCase):

    def test_arbitrary_number_of_sources(self):
        """
        Test that we can define Triple source models.
        """
        n_sources = 16
        parameters = {'t_E': 9}
        for i in range(n_sources):
            parameters['t_0_{0}'.format(i + 1)] = 2. + i
            parameters['u_0_{0}'.format(i + 1)] = 1. / (i + 1)

        model_params = mm.ModelParameters(parameters)
        assert model_params.n_sources == n_sources


class Test1L3SModelErrors(unittest.TestCase):

    def test_incomplete_1L3S_params(self):
        """
        Check error is raised for missing source parameters.
        """
        parameters = {'t_0_1': 0, 'u_0_1': 1,
                      't_0_2': 5, 'rho_2': 0.001,
                      't_0_3': 2, 'u_0_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

        parameters = {'t_0_1': 0, 'u_0_1': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
                      'u_0_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

        parameters = {'t_0_1': 0, 'u_0_1': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
                      'rho_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

    def test_bad_1L3S_params(self):
        """
        Check error is raised if parameter is not associated with a
        specific source.
        """
        parameters = {'t_0': 0, 'u_0_1': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
                      't_0_3': 2, 'u_0_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

        parameters = {'t_0': 0, 'u_0': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
                      't_0_3': 2, 'u_0_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

        parameters = {'t_0': 0, 'u_0_1': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho': 0.001,
                      't_0_3': 2, 'u_0_3': 0.3,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

        parameters = {'t_0': 0, 'u_0_1': 1,
                      't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
                      't_0_3': 2, 'u_0_3': 0.3, 't_star': 0.02,
                      't_E': 9}
        with self.assertRaises(KeyError):
            _ = mm.ModelParameters(parameters)

    def test_1L3S_xallarap(self):
        """
        This test should probably fail, but is it possible to implement this behavior as in OGLE-2015-BLG-1459L?

        https://ui.adsabs.harvard.edu/abs/2018AJ....155..259H/abstract
        """
        parameters = {**xallarap_parameters, 't_0_3': 123.456}
        with self.assertRaises(KeyError):
            mm.ModelParameters(parameters)

        parameters = {**xallarap_parameters, 'u_0_3': 0.456}
        with self.assertRaises(KeyError):
            mm.ModelParameters(parameters)


class TestSetters(unittest.TestCase):

    def setUp(self):
        self.t_0 = 0.0
        self.u_0 = 0.1
        self.t_E = 30
        self.t_eff = self.u_0 * self.t_E
        self.rho = 0.001
        self.t_star = self.rho * self.t_E

        self.q = 0.0001
        self.s = 0.9
        self.alpha = 270.

        self.pi_E_N = 0.5
        self.pi_E_E = -0.3
        self.pi_E = [self.pi_E_N, self.pi_E_E]

        self.dummy_value = 13.

    def test_set_u0(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        params.u_0 = self.dummy_value
        assert params.u_0 == self.dummy_value

    def test_set_u0_error(self):
        params = mm.ModelParameters({'t_0': self.t_0, 't_eff': self.t_eff, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.u_0 = self.dummy_value

    def test_set_t_star_error_1(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 'rho': self.rho})
        with self.assertRaises(AttributeError):
            params.t_star = self.dummy_value

    def test_set_t_star_error_2(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 't_star': self.t_star})
        with self.assertRaises(ValueError):
            params.t_star = -self.dummy_value

    def test_set_t_eff(self):
        params = mm.ModelParameters({'t_0': self.t_0, 't_eff': self.t_eff, 't_E': self.t_E})
        params.t_eff = self.dummy_value
        assert params.t_eff == self.dummy_value

    def test_set_t_eff_error(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.t_eff = self.dummy_value

    def test_set_t_E_error_1(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.t_E = None

    def test_set_t_E_error_2(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.t_E = -self.dummy_value

    def test_set_t_E_error_3(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_eff': self.t_eff})
        with self.assertRaises(AttributeError):
            params.t_E = self.dummy_value

    def test_set_rho(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 'rho': self.rho})
        params.rho = self.dummy_value
        assert params.rho == self.dummy_value

    def test_set_rho_error_1(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 'rho': self.rho})
        with self.assertRaises(ValueError):
            params.rho = -self.dummy_value

    def test_set_rho_error_2(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.rho = self.dummy_value

    def test_set_q(self):
        params = mm.ModelParameters(
            {'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 's': self.s, 'q': self.q, 'alpha': self.alpha})
        params.q = self.dummy_value
        assert params.q == self.dummy_value

    def test_set_q_error(self):
        params = mm.ModelParameters(
            {'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 's': self.s, 'q': self.q, 'alpha': self.alpha})
        with self.assertRaises(ValueError):
            params.q = -self.dummy_value

    def test_set_pi_E(self):
        params = mm.ModelParameters(
            {'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E, 'pi_E_N': self.pi_E_N, 'pi_E_E': self.pi_E_E})
        params.pi_E_N = self.dummy_value
        params.pi_E_E = -self.dummy_value
        assert params.pi_E_N == self.dummy_value
        assert params.pi_E_E == -self.dummy_value

    def test_set_pi_E_error(self):
        params = mm.ModelParameters({'t_0': self.t_0, 'u_0': self.u_0, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.pi_E = [self.dummy_value, -self.dummy_value]


class TestParallax3Sources(unittest.TestCase):

    def setUp(self):
        self.pi_E_N = 0.5
        self.pi_E_E = -0.3
        self.static_params = {'t_0_1': 0, 'u_0_1': 1,
                              't_0_2': 5, 'u_0_2': 0.01,
                              't_0_3': 2, 'u_0_3': 0.3,
                              't_E': 9}
        self.parallax_params = {'pi_E_N': self.pi_E_N, 'pi_E_E': self.pi_E_E}

        self.dummy_value = 13.

    def test_initialized_correctly(self):
        params = mm.ModelParameters({**self.static_params, **self.parallax_params})
        for i in range(3):
            assert params.__getattr__('source_{0}_parameters'.format(i+1)).pi_E_N == self.pi_E_N
            assert params.__getattr__('source_{0}_parameters'.format(i+1)).pi_E_E == self.pi_E_E

    def test_set_pi_E_N(self):
        params = mm.ModelParameters({**self.static_params, 'pi_E_N': self.pi_E_N, 'pi_E_E': self.pi_E_E})
        params.pi_E_N = self.dummy_value
        for i in range(3):
            assert params.__getattr__('source_{0}_parameters'.format(i + 1)).pi_E_N == self.dummy_value
            assert params.__getattr__('source_{0}_parameters'.format(i + 1)).pi_E_E == self.pi_E_E

    def test_set_pi_E_N_error(self):
        params = mm.ModelParameters(self.static_params)
        with self.assertRaises(AttributeError):
            params.pi_E_N = self.dummy_value

    def test_set_pi_E_E(self):
        params = mm.ModelParameters({**self.static_params, 'pi_E_N': self.pi_E_N, 'pi_E_E': self.pi_E_E})
        params.pi_E_E = self.dummy_value
        for i in range(3):
            assert params.__getattr__('source_{0}_parameters'.format(i + 1)).pi_E_N == self.pi_E_N
            assert params.__getattr__('source_{0}_parameters'.format(i + 1)).pi_E_E == self.dummy_value

    def test_set_pi_E_E_error(self):
        params = mm.ModelParameters(self.static_params)
        with self.assertRaises(AttributeError):
            params.pi_E_E = self.dummy_value


class TestSetters2Sources(unittest.TestCase):

    def setUp(self):
        self.t_0_1 = 0.0
        self.u_0_1 = 0.1
        self.t_0_2 = 5.0
        self.u_0_2 = 0.01

        self.t_E = 30
        self.t_eff_1 = self.u_0_1 * self.t_E
        self.t_eff_2 = self.u_0_2 * self.t_E

        self.rho_1 = 0.001
        self.rho_2 = 0.0003
        self.t_star_1 = self.rho_1 * self.t_E
        self.t_star_2 = self.rho_2 * self.t_E

        self.dummy_value = 13.

    def test_set_t_0_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        params.t_0_2 = self.dummy_value
        assert params.t_0_2 == self.dummy_value
        assert params._source_2_parameters.t_0 == self.dummy_value

    def test_set_u_0_1_error(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 't_eff_1': self.t_eff_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.u_0_1 = self.dummy_value

    def test_set_u_0_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        params.u_0_2 = self.dummy_value
        assert params.u_0_2 == self.dummy_value
        assert params._source_2_parameters.u_0 == self.dummy_value

    def test_set_u_0_2_error(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 't_eff_2': self.t_eff_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.u_0_2 = self.dummy_value

    def test_set_rho_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 'rho_1': self.rho_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 'rho_2': self.rho_2, 't_E': self.t_E})
        params.rho_1 = self.dummy_value
        assert params.rho_1 == self.dummy_value
        assert params._source_1_parameters.rho == self.dummy_value

    def test_set_rho_1_error_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.rho_1 = self.dummy_value

    def test_set_rho_1_error_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 'rho_1': self.rho_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 'rho_2': self.rho_2, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.rho_1 = -self.dummy_value

    def test_set_rho_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 'rho_1': self.rho_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 'rho_2': self.rho_2, 't_E': self.t_E})
        params.rho_2 = self.dummy_value
        assert params.rho_2 == self.dummy_value
        assert params._source_2_parameters.rho == self.dummy_value

    def test_set_rho_2_error_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.rho_2 = self.dummy_value

    def test_set_rho_2_error_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 'rho_1': self.rho_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 'rho_2': self.rho_2, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.rho_2 = -self.dummy_value

    def test_set_t_star_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_star_1': self.t_star_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_star_2': self.t_star_2, 't_E': self.t_E})
        params.t_star_1 = self.dummy_value
        assert params.t_star_1 == self.dummy_value
        assert params._source_1_parameters.t_star == self.dummy_value

    def test_set_t_star_1_error_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.t_star_1 = self.dummy_value

    def test_set_t_star_1_error_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_star_1': self.t_star_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_star_2': self.t_star_2, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.t_star_1 = -self.dummy_value

    def test_set_t_star_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_star_1': self.t_star_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_star_2': self.t_star_2, 't_E': self.t_E})
        params.t_star_2 = self.dummy_value
        assert params.t_star_2 == self.dummy_value
        assert params._source_2_parameters.t_star == self.dummy_value

    def test_set_t_star_2_error_1(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_E': self.t_E})
        with self.assertRaises(AttributeError):
            params.t_star_2 = self.dummy_value

    def test_set_t_star_2_error_2(self):
        params = mm.ModelParameters(
            {'t_0_1': self.t_0_1, 'u_0_1': self.u_0_1, 't_star_1': self.t_star_1,
             't_0_2': self.t_0_2, 'u_0_2': self.u_0_2, 't_star_2': self.t_star_2, 't_E': self.t_E})
        with self.assertRaises(ValueError):
            params.t_star_2 = -self.dummy_value


def test_print_3_sources():
    """
    Test if printing of 3-source ModelParameters is correct.
    """
    params = mm.ModelParameters({
        't_0_1': 0, 'u_0_1': 1, 't_E': 9, 't_0_2': 5, 'u_0_2': 0.1, 'rho_2': 0.001,
        't_0_3': 2, 'u_0_3': 0.3, 't_star_3': 0.02})
    expected = ("   t_E (d) \n"
                "    9.0000 \n"
                "  t_0_1 (HJD)     u_0_1 \n"
                "      0.00000  1.000000 \n"
                "  t_0_2 (HJD)     u_0_2   rho_2 \n"
                "      5.00000  0.100000 0.00100 \n"
                "  t_0_3 (HJD)     u_0_3  t_star_3 (d) \n"
                "      2.00000  0.300000      0.020000 ")
    assert params.__repr__() == expected


def _get_times(parameters):
    """prepare a list of epochs based on parameters provided"""
    N = 5
    times = []
    t_0 = parameters.t_0
    t_E = parameters.t_E

    for i in range(N):
        times.append(t_0 - 3 * t_E + i * (6 * t_E / (N - 1)))
    return times


def _get_parameters(orbit):
    """
    Prepare ModelParameters object for testing keplerian orbit calculations
    """
    s = 1.2
    q = 0.123
    u_0 = 0.555
    alpha = 17.5
    t_0 = 2455500
    t_E = 100
    rho = 0.01

    ds_dt = 20.2
    dalpha_dt = -30.3
    s_z = -1.212000
    ds_z_dt = 20
    a = 2

    d = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,
         'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt}

    if orbit == 'circular':
        pass
    elif orbit == 'elliptical':
        d.update({'a_r': a / np.sqrt(s**2+s_z**2), 's_z': s_z})
    else:
        raise KeyError('Requires either circular or elliptical')

    parameters_extra = mm.ModelParameters(d)

    return parameters_extra


def test_separation_for_circular_orbit():
    """
    compares separation to values from VBBinaryLensing v3.6
    """
    parameters = _get_parameters('circular')
    times = _get_times(parameters)

    separation = parameters.get_s(times)

    separation_VBB_circular = [0.57896196, 0.36839104, 1.20000000, 1.66168431, 1.61032808]

    np.testing.assert_almost_equal(separation, separation_VBB_circular)


def test_trajectory_for_circular_orbit():
    """
    compares trajectory to values from VBBinaryLensing v3.6
    """
    parameters = _get_parameters('circular')
    times = _get_times(parameters)

    trajectory = mm.Trajectory(parameters=parameters, times=times)

    x_VBB_circular = [-3.00057308, 1.59071791, 0.16689172, -1.25159957, -2.66309601]
    y_VBB_circular = [-0.55189331, -0.16625747, -0.52931291, -0.99575273, -1.48860493]

    np.testing.assert_almost_equal(trajectory.x, x_VBB_circular)
    np.testing.assert_almost_equal(trajectory.y, y_VBB_circular)


# def test_separation_for_elliptical_orbit():
#     """
#     compares separation to values from VBBinaryLensing v3.6
#     """
#     parameters = _get_parameters('elliptical')
#     times = _get_times(parameters)

#     separation = parameters.get_s(times)

#     separation_VBB_elliptical = [0.64550215, 1.42405587, 1.20000000, 1.34447172, 2.14780495]

#     np.testing.assert_almost_equal(separation, separation_VBB_elliptical)


# def test_trajectory_for_elliptical_orbit():
#     """
#     compares trajectory to values from VBBinaryLensing v3.6
#     """
#     parameters = _get_parameters('elliptical')
#     times = _get_times(parameters)

#     trajectory = mm.Trajectory(parameters=parameters, times=times)

#     x_VBB_elliptical = [-3.03964927, 1.59911138, 0.16689172, 1.23132793, 2.67836603]
#     y_VBB_elliptical = [-0.26183446, -0.02945825, -0.52931291, 1.02071373, 1.46095188]

#     np.testing.assert_almost_equal(trajectory.x, x_VBB_elliptical)
#     np.testing.assert_almost_equal(trajectory.y, y_VBB_elliptical)
