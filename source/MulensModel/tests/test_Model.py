import numpy as np
from numpy.testing import assert_almost_equal as almost
from math import isclose
import unittest
from astropy import units as u

import MulensModel as mm


def test_n_lenses():
    """check n_lenses property"""
    model_1 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.})
    model_2 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.,
                        "s": 1.1234, "q": 0.123, 'alpha': 12.34})
    model_3 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.,
                        "s": 1.1234, "q": 0.123, 'alpha': 12.34,
                        'convergence_K': 0.04, 'shear_G': complex(0.1, -0.05)})
    assert model_1.n_lenses == 1
    assert model_2.n_lenses == 2
    assert model_3.n_lenses == 2


# Point Lens Tests
def test_model_PSPL_1():
    """tests basic evaluation of Paczynski model"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    times = np.array([t_0-2.5*t_E, t_0, t_0+t_E])
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    almost(model.get_magnification(times),
           np.array([1.028720763, 2.10290259, 1.26317278]),
           err_msg="PSPL model returns wrong values")


def test_model_init_1():
    """tests if basic parameters of Model.__init__() are properly passed"""
    t_0 = 5432.10987
    u_0 = 0.001
    t_E = 123.456
    rho = 0.0123
    my_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})
    almost(my_model.parameters.t_0, t_0, err_msg='t_0 not set properly')
    almost(my_model.parameters.u_0, u_0, err_msg='u_0 not set properly')
    almost(my_model.parameters.t_E, t_E, err_msg='t_E not set properly')
    almost(my_model.parameters.rho, rho, err_msg='rho not set properly')


class TestModel(unittest.TestCase):
    def test_negative_t_E(self):
        with self.assertRaises(ValueError):
            my_model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': -100.})


def test_model_parallax_definition():
    """Update parameters in an existing model"""
    model_2 = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.,
                        'pi_E_N': 0.1, 'pi_E_E': 0.2})

    model_2.parameters.pi_E_N = 0.3
    model_2.parameters.pi_E_E = 0.4
    assert model_2.parameters.pi_E_N == 0.3
    assert model_2.parameters.pi_E_E == 0.4

    model_3 = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.,
                        'pi_E': (0.5, 0.6)})
    assert model_3.parameters.pi_E_N == 0.5
    assert model_3.parameters.pi_E_E == 0.6

    model_4 = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.,
                        'pi_E_N': 0.7, 'pi_E_E': 0.8})
    assert model_4.parameters.pi_E_N == 0.7
    assert model_4.parameters.pi_E_E == 0.8


def test_coords_transformation():
    """
    this was tested using http://ned.ipac.caltech.edu/forms/calculator.html
    """
    coords = "17:54:32.1 -30:12:34.0"
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.}, coords=coords)

    almost(model.coords.galactic_l.value, 359.90100049-360., decimal=4)
    almost(model.coords.galactic_b.value, -2.31694073, decimal=3)

    almost(model.coords.ecliptic_lon.value, 268.81102051, decimal=1)
    almost(model.coords.ecliptic_lat.value, -6.77579203, decimal=2)


def test_init_parameters():
    """are parameters properly passed between Model and ModelParameters?"""
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    model = mm.Model(parameters=params)
    almost(model.parameters.t_0, t_0)
    almost(model.parameters.u_0, u_0)
    almost(model.parameters.t_E, t_E.value)


def test_limb_darkening():
    """check if limb_darkening coeffs are properly passed and converted"""
    gamma = 0.4555
    u = 3. * gamma / (2. + gamma)

    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'rho': 0.001})
    model.set_limb_coeff_gamma('I', gamma)

    almost(model.get_limb_coeff_gamma('I'), gamma)
    almost(model.get_limb_coeff_u('I'), u)


def test_t_E():
    """make sure t_E can be accessed properly"""
    t_0 = 2460000.
    u_0 = 0.1
    t_E = 12.3456
    t_star = 0.01234
    params = dict(t_0=t_0, u_0=u_0, t_E=t_E)
    model_1 = mm.Model(params)
    params['t_star'] = t_star
    model_2 = mm.Model(params)

    almost(model_1.parameters.t_E, t_E)
    almost(model_2.parameters.t_E, t_E)


# Binary Lens tests
# Binary lens parameters:
alpha = 229.58 * u.deg
s = 1.3500
q = 0.00578
# Other parameters
t_E = 62.63 * u.day
rho = 0.01

# Adjust t_0, u_0 from primary to CM
shift_x = -s * q / (1. + q)
delta_t_0 = -t_E.value * shift_x * np.cos(alpha).value
delta_u_0 = -shift_x * np.sin(alpha).value
t_0 = 2456141.593 + delta_t_0
u_0 = 0.5425 + delta_u_0


def test_BLPS_01():
    """simple binary lens with point source"""
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q})

    model = mm.Model(parameters=params)
    t = np.array([2456112.5])
    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    magnification = model.get_magnification(data.time[0])
    almost(magnification, 4.691830781584699)
# This value comes from early version of this code.
# almost(m, 4.710563917)
# This value comes from Andy's getbinp().


def test_BLPS_shear_active():
    """
    simple binary lens with point source and external convergence and shear
    """
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'convergence_K': 0.08, 'shear_G': complex(0.1, -0.1)})

    model = mm.Model(parameters=params)
    t = np.array([2456112.5])
    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    magnification = model.get_magnification(data.time[0])
    assert not isclose(magnification, 4.691830781584699, abs_tol=1e-2)


def test_BLPS_shear():
    """
    simple binary lens with point source and external convergence and shear
    """
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'convergence_K': 0.0, 'shear_G': complex(0, 0)})

    model = mm.Model(parameters=params)
    t = np.array([2456112.5])
    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    magnification = model.get_magnification(data.time[0])
    almost(magnification, 4.691830781584699)


def test_BLPS_02():
    """
    simple binary lens with extended source and different methods to
    evaluate magnification
    """

    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'rho': rho})
    model = mm.Model(parameters=params)

    t = np.array([6112.5, 6113., 6114., 6115., 6116., 6117., 6118., 6119])
    t += 2450000.
    methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole', 2456116.5,
               'VBBL', 2456117.5]
    model.set_magnification_methods(methods)

    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    result = model.get_magnification(data.time)

    expected = np.array([4.69183078, 2.87659723, 1.83733975, 1.63865704,
                         1.61038135, 1.63603122, 1.69045492, 1.77012807])
    almost(result, expected)

    # Possibly, this test should be re-created in test_FitData.py
    # Below we test passing the limb coeff to VBBL function.
    # data.bandpass = 'I'
    model.set_limb_coeff_u('I', 10.)
    # This is an absurd value but I needed something quick.
    result = model.get_magnification(
        data.time, gamma=model.get_limb_coeff_gamma('I'))
    almost(result[5], 1.6366862)
    result_2 = model.get_magnification(data.time, bandpass='I')
    almost(result, result_2)


def test_BLPS_02_AC():
    """
    simple binary lens with extended source and different methods to evaluate
    magnification - version with adaptivecontouring
    """

    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'rho': rho})
    model = mm.Model(parameters=params)

    t = np.array([6112.5, 6113., 6114., 6115., 6116., 6117., 6118., 6119])
    t += 2450000.
    ac_name = 'Adaptive_Contouring'
    methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole', 2456116.5,
               ac_name, 2456117.5]
    accuracy_1 = {'accuracy': 0.04}
    accuracy_2 = {'accuracy': 0.01, 'ld_accuracy': 0.00001}
    model.set_magnification_methods(methods)
    model.set_magnification_methods_parameters({ac_name: accuracy_1})

    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    result = model.get_magnification(data.time)

    expected = np.array([4.69183078, 2.87659723, 1.83733975, 1.63865704,
                         1.61038135, 1.63603122, 1.69045492, 1.77012807])
    almost(result, expected, decimal=3)

    # Below we test passing the limb coeff to VBBL function.
    # data.bandpass = 'I'
    model.set_limb_coeff_u('I', 10.)
    # This is an absurd value but I needed something quick.
    model.set_magnification_methods_parameters({ac_name: accuracy_2})
    # result = model.data_magnification[0]
    result = model.get_magnification(
        data.time, gamma=model.get_limb_coeff_gamma('I'))
    almost(result[5], 1.6366862, decimal=3)


def test_methods_parameters():
    """
    make sure additional parameters are properly passed to very inner functions
    """
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'rho': rho})
    model = mm.Model(parameters=params)

    t = np.array([2456117.])
    methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole', 2456116.5,
               'VBBL', 2456117.5]
    model.set_magnification_methods(methods)

    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    result_1 = model.get_magnification(data.time)

    vbbl_options = {'accuracy': 0.1}
    methods_parameters = {'VBBL': vbbl_options}
    model.set_magnification_methods_parameters(methods_parameters)
    result_2 = model.get_magnification(data.time)

    vbbl_options = {'accuracy': 1.e-5}
    methods_parameters = {'VBBL': vbbl_options}
    model.set_magnification_methods_parameters(methods_parameters)
    result_3 = model.get_magnification(data.time)

    assert result_1[0] != result_2[0]
    assert result_1[0] != result_3[0]
    assert result_2[0] != result_3[0]


def test_caustic_for_orbital_motion():
    """
    check if caustics calculated for different epochs in orbital motion model
    are as expected
    """
    q = 0.1
    s = 1.3
    params = {'t_0': 100., 'u_0': 0.1, 't_E': 10., 'q': q,
              's': s, 'ds_dt': 0.5, 'alpha': 0., 'dalpha_dt': 0.}
    model = mm.Model(parameters=params)

    model.update_caustics()
    almost(model.caustics.get_caustics(),
           mm.Caustics(q=q, s=s).get_caustics())

    model.update_caustics(100.+365.25/2)
    almost(model.caustics.get_caustics(),
           mm.Caustics(q=q, s=1.55).get_caustics())


def test_magnifications_for_orbital_motion():
    """
    make sure that orbital motion parameters are properly passed to
    magnification methods calculations
    """
    dict_static = {'t_0': 100., 'u_0': 0.1, 't_E': 100., 'q': 0.99,
                   's': 1.1, 'alpha': 10.}
    dict_motion = dict_static.copy()
    dict_motion.update({'ds_dt': -2, 'dalpha_dt': -300.})
    static = mm.Model(dict_static)
    motion = mm.Model(dict_motion)

    t_1 = 100.
    almost(static.get_magnification(t_1), motion.get_magnification(t_1))

    t_2 = 130.
    static.parameters.s = 0.93572895
    static.parameters.alpha = 345.359342916
    almost(static.get_magnification(t_2), motion.get_magnification(t_2))


def test_model_binary_and_finite_sources():
    """
    test if model magnification calculation for binary source works with
    finite sources (both rho and t_star given)
    """
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.005, 'rho_1': 0.001,
        't_0_2': 5100., 'u_0_2': 0.0003, 't_star_2': 0.03, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    (t1, t2) = (4999.95, 5000.05)
    (t3, t4) = (5099.95, 5100.05)
    finite = 'finite_source_uniform_Gould94'
    model.set_magnification_methods(
        [t1, finite, t2, 'point_source', t3, finite, t4])
    model_1.set_magnification_methods([t1, finite, t2])
    model_2.set_magnification_methods([t3, finite, t4])

    (f_s_1, f_s_2, f_b) = (100., 300., 50.)
    time = np.linspace(4900., 5200., 4200)
    mag_1 = model_1.get_magnification(time)
    mag_2 = model_2.get_magnification(time)

    # test: model.set_source_flux_ratio(f_s_2/f_s_1)
    fitted = model.get_magnification(time, source_flux_ratio=f_s_2 / f_s_1)
    expected = (mag_1 * f_s_1 + mag_2 * f_s_2) / (f_s_1 + f_s_2)
    almost(fitted, expected)

    # test separate=True option:
    (mag_1_, mag_2_) = model.get_magnification(time, separate=True)
    almost(mag_1, mag_1_)
    almost(mag_2, mag_2_)


def test_binary_source_and_fluxes_for_bands():
    """
    Test if setting different flux ratios for different bands in binary
    source models works properly. The argument flux_ratio_constraint
    is set as string.
    """
    model = mm.Model({'t_0_1': 5000., 'u_0_1': 0.05,
                      't_0_2': 5100., 'u_0_2': 0.003, 't_E': 30.})

    times_I = np.linspace(4900., 5200., 3000)
    times_V = np.linspace(4800., 5300., 250)
    (f_s_1_I, f_s_2_I, f_b_I) = (10., 20., 3.)
    (f_s_1_V, f_s_2_V, f_b_V) = (15., 5., 30.)
    q_f_I = f_s_2_I / f_s_1_I
    q_f_V = f_s_2_V / f_s_1_V
    (mag_1_I, mag_2_I) = model.get_magnification(times_I, separate=True)
    (mag_1_V, mag_2_V) = model.get_magnification(times_V, separate=True)
    effective_mag_I = (mag_1_I + mag_2_I * q_f_I) / (1. + q_f_I)
    effective_mag_V = (mag_1_V + mag_2_V * q_f_V) / (1. + q_f_V)
    flux_I = mag_1_I * f_s_1_I + mag_2_I * f_s_2_I + f_b_I
    flux_V = mag_1_V * f_s_1_V + mag_2_V * f_s_2_V + f_b_V

    # model.set_source_flux_ratio_for_band('I', q_f_I)
    # model.set_source_flux_ratio_for_band('V', q_f_V)

    # Test Model.get_magnification()
    result_I = model.get_magnification(times_I, source_flux_ratio=q_f_I)
    result_V = model.get_magnification(times_V, source_flux_ratio=q_f_V)
    almost(result_I, effective_mag_I)
    almost(result_V, effective_mag_V)


def test_separate_method_for_each_source():
    """
    Checks if setting separate magnification method for each source in
    binary source models works properly.
    """
    model = mm.Model({'t_0_1': 5000., 'u_0_1': 0.01, 'rho_1': 0.005,
                      't_0_2': 5100., 'u_0_2': 0.001, 'rho_2': 0.005,
                      't_E': 1000.})

    model.set_magnification_methods(
        [4999., 'finite_source_uniform_Gould94', 5101.], source=1)
    # In order not to get "no FS method" warning:
    model.set_magnification_methods(
        [0., 'finite_source_uniform_Gould94', 1.], source=2)
    out = model.get_magnification(5000., separate=True)
    almost([out[0][0], out[1][0]], [103.46704167, 10.03696291])

    model.set_magnification_methods(
        [4999., 'finite_source_uniform_Gould94', 5001.], source=1)
    model.set_magnification_methods(
        [5099., 'finite_source_uniform_Gould94', 5101.], source=2)
    out = model.get_magnification(5100., separate=True)
    almost([out[0][0], out[1][0]], [9.98801936, 395.96963727])


def test_get_lc():
    """
    Test if Model.get_lc() works properly; we test on binary source model
    without finite source effect.
    """
    model = mm.Model({'t_0_1': 5000., 'u_0_1': 1.,
                      't_0_2': 5100., 'u_0_2': 0.1,
                      't_E': 100.})
    out = model.get_lc(5050., source_flux=[1., 2.], blend_flux=3.)
    almost(out, 19.668370500043526)


# Tests to Add:
#
# test get_trajectory:
#   straight-up trajectory
#   case with annual parallax (check test_Model_Parallax.py)
#   case with satellite parallax (check test_Model_Parallax.py)
#   coords is propagating correctly (check test_Model_Parallax.py)
#
# test set_times:
#   keywords to test:
#     t_range=None, t_start=None, t_stop=None, dt=None, n_epochs=1000
#
# test set_default_magnification_method:
#   change from default value
#
# test get_satellite_coords: (check test_Model_Parallax.py)
#   returns None if no ephemerides file set
#   other condidtions probably covered by other unit tests
#
# test _magnification_2_sources: check instances of q_flux being specified vs.
# separate. Specifically worried that calls to magnification from other parts
# of the code work as expected.
#
# properties: parallax, caustics, parameters, n_lenses, n_source, is_static,
# coords, bandpasses,
