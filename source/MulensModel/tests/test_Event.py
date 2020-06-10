import sys
from os.path import join
import unittest
import numpy as np
from astropy import units as u

import MulensModel as mm


dir_ = join(mm.DATA_PATH, "photometry_files")
SAMPLE_FILE_01 = join(dir_, "OB08092", "phot_ob08092_O4.dat")
SAMPLE_FILE_02 = join(dir_, "OB140939", "ob140939_OGLE.dat")
SAMPLE_FILE_03 = join(dir_, "OB03235", "OB03235_OGLE.tbl.txt")
SAMPLE_FILE_04 = join(dir_, "OB03235", "OB03235_MOA.tbl.txt")
SAMPLE_FILE_310_01 = join(dir_, "MB08310", "MOA_0300089_PLC_007.tbl")
SAMPLE_FILE_310_02 = join(dir_, "MB08310", "Bron_0300089_PLC_002.tbl")
SAMPLE_FILE_310_03 = join(dir_, "MB08310", "CTIO_I_0300089_PLC_005.tbl")

def test_model_event_coords():
    """
    Check if coordinates are properly passed from Model to Event.
    """
    coords = "18:12:34.56 -23:45:55.55"
    model = mm.Model({'t_0': 0, 'u_0': .5, 't_E': 10.}, coords=coords)
    event = mm.Event(model=model)
    np.testing.assert_almost_equal(event.coords.ra.value, 273.144)
    np.testing.assert_almost_equal(event.coords.dec.value, -23.765430555555554)


def test_event_get_chi2_1():
    """basic unit test on ob08092 OGLE-IV data"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data])
    ev.model = mod
    ev.datasets = [data]

    # Make sure Event.fit is defined (can be None):
    assert ev.fit is None or isinstance(ev.fit, mm.Fit)

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 427.20382, decimal=4,
                                   err_msg='problem in resulting chi2')

    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')

    ev.sum_function = 'numpy.sum'
    np.testing.assert_almost_equal(
        ev.get_chi2(), 427.20382, decimal=4, err_msg='problem with numpy.sum')


def test_event_get_chi2_2():
    """
    Basic unit test on ob08092 OGLE-IV data. Same as above but with
    the data input twice (to test behavior for multiple datasets);
    also test if Event.get_chi2_for_dataset() gives correct answers.
    """
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    answer = 427.20382201

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data, data])
    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 2.*answer, decimal=4,
                                   err_msg='problem in resulting chi2')

    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 2.*459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')

    chi2_2 = ev.get_chi2_for_dataset(0)
    np.testing.assert_almost_equal(chi2_2, answer)

    chi2_3 = ev.get_chi2_for_dataset(1)
    np.testing.assert_almost_equal(chi2_3, answer)


def test_event_get_chi2_3():
    """test on ob08092 OGLE-IV data - MulensData.good & MulensData.bad"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    bad = np.zeros(383, dtype='bool')
    bad[300:350] = True
    data_1 = mm.MulensData(file_name=SAMPLE_FILE_01, bad=bad)
    data_2 = mm.MulensData(file_name=SAMPLE_FILE_01, good=~bad)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data_1])
    ev.model = mod
    ev.datasets = [data_1]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 343.46567, decimal=4,
                                   err_msg='problem in resulting chi2')

    mod.set_datasets([data_2])
    ev.model = mod
    ev.datasets = [data_2]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 343.46567, decimal=4,
                                   err_msg='problem in resulting chi2')


def test_event_get_chi2_double_source_simple():
    """
    basic test on ob08092 OGLE-IV data
    """
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data, data])

    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()

    assert isinstance(chi2, float), 'wrong type of chi2'
    message = 'problem in resulting chi2 for 2 exactly the same datasets'
    np.testing.assert_almost_equal(
        chi2, 854.407644, decimal=4, err_msg=message)


def test_event_get_chi2_3():
    """
    Test: If I change the model parameters, the chi2 should change.
    """
    t_0 = 5380.
    u_0 = 0.5
    t_E = 18.
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

    # Generate fake data
    times = np.arange(5320, 5420.)
    f_source = 0.1
    f_blend = 0.5
    mod_fluxes = f_source * model.magnification(times) + f_blend
    I_mag = mm.Utils.get_mag_from_flux(mod_fluxes)
    errors = 0.01 * np.ones(times.shape)
    data = mm.MulensData(data_list=[times, I_mag, errors])

    # Generate event and fit
    event = mm.Event(datasets=data, model=model)

    orig_chi2 = event.get_chi2()

    # Change the model
    event.model.parameters.t_0 = 5000.

    assert event.get_chi2() != orig_chi2


def test_event_get_chi2_4():
    """test if best chi2 is remembered correctly"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E*u.day}
    mod = mm.Model(params)
    mod.set_datasets([data])
    ev.model = mod
    ev.datasets = [data]

    chi2_1 = ev.get_chi2()  # This is the best model.

    ev.model.parameters.parameters['t_0'] += 1.
    ev.model.parameters.parameters['u_0'] += 0.1
    ev.model.parameters.parameters['t_E'] += 1. * u.day
    chi2_2 = ev.get_chi2()

    assert chi2_2 > chi2_1
    assert ev.best_chi2 == chi2_1
    assert ev.best_chi2_parameters == params


def test_event_get_chi2_5():
    """
    De facto checks how information is passed between Event, Fit, and Model.
    Also checks simple relations between different Event.get_chi2().
    """
    kwargs = {'comments': ["\\", "|"]}
    data_1 = mm.MulensData(file_name=SAMPLE_FILE_03, **kwargs)
    data_2 = mm.MulensData(file_name=SAMPLE_FILE_04, phot_fmt='flux', **kwargs)
    params = {'t_0': 2452848.06, 'u_0': 0.1317, 't_E': 61.5}

    model = mm.Model(params)

    event_1 = mm.Event(data_1, model)
    event_2 = mm.Event(data_2, model)
    event_1_2 = mm.Event([data_1, data_2], model)
    event_2_1 = mm.Event([data_2, data_1], model)

    chi2_1 = event_1.get_chi2()
    chi2_2 = event_2.get_chi2()
    # The line below was failing when the test was written:
    chi2_1_2 = event_1_2.get_chi2()
    chi2_2_1 = event_2_1.get_chi2()

    np.testing.assert_almost_equal(chi2_1_2, chi2_1 + chi2_2)
    np.testing.assert_almost_equal(chi2_1_2, chi2_2_1)


class TestEvent(unittest.TestCase):
    def test_event_init_1(self):
        with self.assertRaises(TypeError):
            ev = mm.Event(model=3.14)

    def test_event_init_2(self):
        with self.assertRaises(TypeError):
            ev = mm.Event(datasets='some_string')


def test_event_chi2_gradient():
    """test calculation of chi2 gradient"""
    # fs = 11.0415734, fb = 0.0
    parameters_1 = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87}
    params_1 = ['t_0', 'u_0', 't_E']
    gradient_1 = {'t_0': 236.206598, 'u_0': 101940.249,
                  't_E': -1006.88678}
    test_1 = (parameters_1, params_1, gradient_1)

    parameters_2 = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
                    'pi_E_N': -0.248, 'pi_E_E': 0.234}
    # This model also used fluxes given above.
    params_2 = ['t_0', 'u_0', 't_E', 'pi_E_N', 'pi_E_E', 'f_source', 'f_blend']
    gradient_2 = {'t_0': 568.781786, 'u_0': 65235.3513, 't_E': -491.782005,
                  'pi_E_N': -187878.357, 'pi_E_E': 129162.927,
                  'f_source': -83124.5869, 'f_blend': -78653.242}
    test_2 = (parameters_2, params_2, gradient_2)
    # We're not applying the test above, yet. See 'for' loop below.

    data = mm.MulensData(file_name=SAMPLE_FILE_02)
    kwargs = {'datasets': [data], 'coords': '17:47:12.25 -21:22:58.7'}

    for test in [test_1]:  # , test_2]:
        (parameters, params, gradient) = test
        event = mm.Event(model=mm.Model(parameters), **kwargs)
        result = event.chi2_gradient(params, fit_blending=False)

        reference = np.array([gradient[key] for key in params])
        np.testing.assert_almost_equal(reference/result, 1., decimal=1)

def _test_event_chi2_gradient_rho():
    """test calculation of chi2 gradient including finite source effects"""
    # MB08310 is the example
    # Set datasets
    datasets = [mm.MulensData(file_name=SAMPLE_FILE_310_01, bandpass='R'),
                mm.MulensData(file_name=SAMPLE_FILE_310_02, bandpass='U'),
                mm.MulensData(file_name=SAMPLE_FILE_310_03, bandpass='I')]

    # Setup Model
    gamma_I, gamma_V = (0.44, 0.72)
    t_0 = 2454656.39975
    u_0 = 0.00300
    t_E = 11.14
    rho = 0.00492549
    t_star = rho * t_E
    parameters = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho}
    model = mm.Model(parameters)
    method = 'finite_source_LD_Yoo04'
    model.set_magnification_methods(
        [t_0 - 2. * t_star, method, t_0 + 2. * t_star])
    model.set_limb_coeff_gamma('R', (gamma_V + gamma_I) / 2.)
    model.set_limb_coeff_gamma('U', (gamma_V + gamma_I) / 2.)
    model.set_limb_coeff_gamma('I', gamma_I)

    # Set expected values
    # JCY - see sandbox/rho_gradient on pink laptop
    params = parameters.keys
    gradient = {'t_0': 1283513.3068849628, 'u_0': 20492801.742886964,
                't_E': -9573.3589902395161, 'rho': -1503911.2409404013}
    reference = np.array([gradient[key] for key in params])

    # Create event and run test
    event = mm.Event(model=model, datasets=datasets)
    #result = event.chi2_gradient(params, fit_blending=False)

    #print(result)
    #print(reference)
    #np.testing.assert_almost_equal(reference / result, 1., decimal=2)


def test_get_ref_fluxes():
    """Test Event.get_ref_fluxes()"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    event = mm.Event(data, model)

    (f_s_1, f_b_1) = event.get_ref_fluxes()
    (f_s_2, f_b_2) = event.get_ref_fluxes(fit_blending=False)
    (f_s_3, f_b_3) = event.get_ref_fluxes(fit_blending=True)

    assert f_b_2 == 0.
    assert f_s_1 == f_s_3
    assert f_b_1 == f_b_3
    np.testing.assert_almost_equal((f_s_1 + f_b_1)/f_s_2, 1., decimal=3)
    # Table 1 of Poleski et al. 2014:
    np.testing.assert_almost_equal(f_b_1 / f_s_1, 0.016, decimal=3)


def test_event_chi2_binary_source():
    """simple test if chi2 calculation for binary source works fine"""
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.05,
        't_0_2': 5100., 'u_0_2': 0.15, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    # prepare fake data:
    time = np.linspace(4900., 5200., 600)
    mag_1 = model_1.magnification(time)
    mag_2 = model_2.magnification(time)
    flux = 100. * mag_1 + 300. * mag_2 + 50.
    data = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')

    # Calculate chi^2:
    event = mm.Event([data], model)
    np.testing.assert_almost_equal(event.get_chi2(), 0.)
    # Make sure Model.set_source_flux_ratio() is taken into account.
    model.set_source_flux_ratio(1.)
    np.testing.assert_almost_equal(model.magnification(time), (mag_1+mag_2)/2.)
    assert event.get_chi2() > 1.
    model.set_source_flux_ratio(3.)
    np.testing.assert_almost_equal(event.get_chi2(), 0.)


def test_event_chi2_binary_source_2datasets():
    """
    simple test if chi2 calculation for binary source
    works fine for 2 datasets
    """
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.05,
        't_0_2': 5100., 'u_0_2': 0.15, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    # prepare fake data:
    time = np.linspace(4900., 5200., 600)
    mag_1 = model_1.magnification(time)
    mag_2 = model_2.magnification(time)
    flux = 100. * mag_1 + 300. * mag_2 + 50.
    data_1 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')
    flux = 20. * mag_1 + 30. * mag_2 + 50.
    data_2 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')

    # Calculate chi^2:
    event = mm.Event([data_1, data_2], model)
    np.testing.assert_almost_equal(event.get_chi2(), 0.)
    np.testing.assert_almost_equal(event.get_chi2_for_dataset(0), 0.)
    np.testing.assert_almost_equal(event.get_chi2_for_dataset(1), 0.)
    # Make sure that changing parameters changes chi2:
    model.parameters.t_E = 100.
    assert event.get_chi2() > 1., 'wrong chi2'
    model.parameters.t_E = 25.
    np.testing.assert_almost_equal(event.get_chi2(), 0.)
    model.parameters.t_0_1 = 5010.
    assert event.get_chi2() > 1., 'wrong chi2'
    model.parameters.t_0_1 = 5000.

    # Test combination of Model.set_source_flux_ratio_for_band() and
    # Event.get_chi2_for_dataset().
    data_1.bandpass = 'some'
    event.model.set_source_flux_ratio_for_band('some', 3.)
    np.testing.assert_almost_equal(event.get_chi2_for_dataset(0), 0.)
