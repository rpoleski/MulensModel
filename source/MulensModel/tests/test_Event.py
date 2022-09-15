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


# ---------
# Functions used by unit tests
def generate_binary_source_models():
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.05,
        't_0_2': 5100., 'u_0_2': 0.15, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    return (model, model_1, model_2)


def generate_binary_source_datasets(model_1, model_2):
    """prepare fake data for binary source model"""
    time = np.linspace(4900., 5200., 600)
    mag_1 = model_1.get_magnification(time)
    mag_2 = model_2.get_magnification(time)
    flux = 100. * mag_1 + 300. * mag_2 + 50.
    data_1 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')
    flux = 20. * mag_1 + 30. * mag_2 + 50.
    data_2 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')

    return (data_1, data_2)


# ----------
# Event Coordinates Tests
def test_model_event_coords():
    """
    Check if coordinates are properly passed from Model to Event.
    """
    coords = "18:12:34.56 -23:45:55.55"
    model = mm.Model({'t_0': 0, 'u_0': .5, 't_E': 10.}, coords=coords)
    event = mm.Event(model=model)
    np.testing.assert_almost_equal(event.coords.ra.value, 273.144)
    np.testing.assert_almost_equal(event.coords.dec.value, -23.765430555555554)


# ----------
# Event.get_chi2() Tests - Basic
def test_event_get_chi2_1():
    """basic unit test on ob08092 OGLE-IV data"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    ev.model = mod
    ev.datasets = [data]

    # Make sure Event.fit is defined (can be None):
    # assert ev.fit is None or isinstance(ev.fit, mm.Fit)

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 427.20382, decimal=4,
                                   err_msg='problem in resulting chi2')

    ev.sum_function = 'numpy.sum'
    np.testing.assert_almost_equal(
        ev.get_chi2(), 427.20382, decimal=4, err_msg='problem with numpy.sum')

    # New (v2.0) method of fixing the blending
    ev.fix_blend_flux = {}
    chi2_no_blend = ev.get_chi2()
    np.testing.assert_almost_equal(chi2_no_blend, 427.20382, decimal=4)

    fix_blend_flux = {}
    for dataset in ev.datasets:
        fix_blend_flux[dataset] = 0.

    ev.fix_blend_flux = fix_blend_flux
    np.testing.assert_almost_equal(ev.get_chi2(), 459.09826, decimal=4)
    assert(ev.blend_fluxes == [0.])


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

    data_1 = mm.MulensData(file_name=SAMPLE_FILE_01)
    data_2 = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    ev.model = mod
    ev.datasets = [data_1, data_2]

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 2.*answer, decimal=4,
                                   err_msg='problem in resulting chi2')

    chi2_2 = ev.get_chi2_for_dataset(0)
    np.testing.assert_almost_equal(chi2_2, answer)

    chi2_3 = ev.get_chi2_for_dataset(1)
    np.testing.assert_almost_equal(chi2_3, answer)

    # New method of fixing the blending
    # Both datasets have zero blending
    ev.fix_blend_flux = {data_1: 0., data_2: 0.}
    chi2_no_blend = ev.get_chi2()
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 2.*459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')
    assert (ev.blend_fluxes == [0., 0.])


def test_event_get_chi2_3():
    """test on ob08092 OGLE-IV data - MulensData.good & MulensData.bad"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    bad = np.zeros(383, dtype='bool')
    bad[300:350] = True

    data_1 = mm.MulensData(file_name=SAMPLE_FILE_01, bad=bad)
    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    ev.model = mod
    ev.datasets = [data_1]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 342.60457, decimal=2,
                                   err_msg='problem in resulting chi2')

    data_2 = mm.MulensData(file_name=SAMPLE_FILE_01, good=~bad)
    ev.model = mod
    ev.datasets = [data_2]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 342.60457, decimal=2,
                                   err_msg='problem in resulting chi2')


def test_event_get_chi2_4():
    """test if best chi2 is remembered correctly"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E * u.day}
    mod = mm.Model(params)
    ev.model = mod
    ev.datasets = [data]

    chi2_1 = ev.get_chi2()  # This is the best model.

    ev.model.parameters.parameters['t_0'] += 1.
    ev.model.parameters.parameters['u_0'] += 0.1
    ev.model.parameters.parameters['t_E'] += 1. * u.day
    chi2_2 = ev.get_chi2()

    assert chi2_2 > chi2_1


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
    chi2_1_2 = event_1_2.get_chi2()
    chi2_2_1 = event_2_1.get_chi2()

    np.testing.assert_almost_equal(chi2_1_2, chi2_1 + chi2_2)
    np.testing.assert_almost_equal(chi2_1_2, chi2_2_1)


def test_event_get_chi2_6():
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
    mod_fluxes = f_source * model.get_magnification(times) + f_blend
    I_mag = mm.Utils.get_mag_from_flux(mod_fluxes)
    errors = 0.01 * np.ones(times.shape)
    data = mm.MulensData(data_list=[times, I_mag, errors])

    # Generate event and fit
    event = mm.Event(datasets=data, model=model)
    orig_chi2 = event.get_chi2()

    # Change the model
    event.model.parameters.t_0 = 5000.

    assert event.get_chi2() != orig_chi2


# ----------
# Event.get_chi2() Tests - Alternate Cases
def test_chi2_vs_get_chi2():
    """get_chi2 always updates after something changes in the model, chi2 does
    not."""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002

    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    ev = mm.Event()
    mod = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    ev.model = mod
    ev.datasets = [data]

    # New (v2.0) method of fixing the blending
    ev.fix_blend_flux = {}
    np.testing.assert_almost_equal(ev.get_chi2(), 427.20382, decimal=4)
    np.testing.assert_almost_equal(ev.chi2, 427.20382, decimal=4)

    fix_blend_flux = {}
    for dataset in ev.datasets:
        fix_blend_flux[dataset] = 0.

    ev.fix_blend_flux = fix_blend_flux
    np.testing.assert_almost_equal(ev.chi2, 427.20382, decimal=4)
    np.testing.assert_almost_equal(ev.get_chi2(), 459.09826, decimal=4)
    np.testing.assert_almost_equal(ev.chi2, 459.09826, decimal=4)


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

    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()

    assert isinstance(chi2, float), 'wrong type of chi2'
    message = 'problem in resulting chi2 for 2 exactly the same datasets'
    np.testing.assert_almost_equal(
        chi2, 854.407644, decimal=4, err_msg=message)


def test_event_chi2_binary_source():
    """simple test if chi2 calculation for binary source works fine"""
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.05,
        't_0_2': 5100., 'u_0_2': 0.15, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    # prepare fake data:
    time = np.linspace(4900., 5200., 600)
    mag_1 = model_1.get_magnification(time)
    mag_2 = model_2.get_magnification(time)
    flux = 100. * mag_1 + 300. * mag_2 + 50.
    data = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')

    # Calculate chi^2:
    event = mm.Event([data], model)
    np.testing.assert_almost_equal(event.get_chi2(), 0.)

    # Make sure fix_source_flux_ratio is taken into account.
    event.fix_source_flux_ratio = {data: 1.}
    assert event.get_chi2() > 1.
    event.fix_source_flux_ratio = {data: 3.}
    np.testing.assert_almost_equal(event.get_chi2(), 0.)


def test_event_chi2_binary_source_2datasets():
    """
    simple test if chi2 calculation for binary source
    works fine for 2 datasets
    """
    (model, model_1, model_2) = generate_binary_source_models()
    (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)

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
    event.fix_source_flux_ratio['some'] = 3.
    event.fit_fluxes()
    np.testing.assert_almost_equal(event.get_chi2_for_dataset(0), 0.)


# --------
# Error message tests
class TestEvent(unittest.TestCase):
    def test_event_init_1(self):
        with self.assertRaises(TypeError):
            _ = mm.Event(model=3.14)

    def test_event_init_2(self):
        with self.assertRaises(TypeError):
            _ = mm.Event(datasets='some_string')


# ----------
# Chi2 Gradient Tests
class Chi2GradientTest():

    def __init__(self, parameters=None, grad_params=None, gradient=None):
        self.parameters = parameters
        self.grad_params = grad_params
        self.gradient = gradient


chi2_gradient_test_1 = Chi2GradientTest(
    parameters={'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87},
    grad_params=['t_0', 'u_0', 't_E'],
    gradient={'t_0': 236.206598, 'u_0': 101940.249, 't_E': -1006.88678})

# Not used:
# f_source and f_blend cannot be gradient parameters in MulensModel, but
# this test could be moved to sfit_minimizer, which is under development by JCY.
chi2_gradient_test_2 = Chi2GradientTest(
    parameters={'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
                'pi_E_N': -0.248, 'pi_E_E': 0.234},
    grad_params=['t_0', 'u_0', 't_E', 'pi_E_N', 'pi_E_E', 'f_source',
                 'f_blend'],
    gradient={'t_0': 568.781786, 'u_0': 65235.3513, 't_E': -491.782005,
              'pi_E_N': -187878.357, 'pi_E_E': 129162.927,
              'f_source': -83124.5869, 'f_blend': -78653.242})


def test_event_chi2_gradient():
    """test calculation of chi2 gradient"""

    data = mm.MulensData(file_name=SAMPLE_FILE_02)
    kwargs = {'datasets': [data], 'coords': '17:47:12.25 -21:22:58.7'}

    # New (v2.0) method for fixing blending
    for test in [chi2_gradient_test_1]:
        event = mm.Event(
            model=mm.Model(test.parameters), fix_blend_flux={data: 0.},
            **kwargs)
        result = event.get_chi2_gradient(test.grad_params)
        reference = np.array([test.gradient[key] for key in test.grad_params])
        np.testing.assert_almost_equal(reference/result, 1., decimal=4)

        result = event.chi2_gradient
        np.testing.assert_almost_equal(reference/result, 1., decimal=4)


class TestGradient(unittest.TestCase):

    def test_gradient_init_1(self):
        """
        Test that fit_fluxes/update must be run before the gradient is accessed
        """
        data = mm.MulensData(file_name=SAMPLE_FILE_02)
        event = mm.Event(
            datasets=[data],
            model=mm.Model(chi2_gradient_test_1.parameters),
            fix_blend_flux={data: 0.})
        with self.assertRaises(TypeError):
            event.chi2_gradient()


def test_chi2_gradient():
    """
    test that fit_fluxes/update must be run in order to update the gradient
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_02)
    event = mm.Event(
        datasets=[data],
        model=mm.Model(chi2_gradient_test_1.parameters),
        fix_blend_flux={data: 0.})
    result = event.get_chi2_gradient(chi2_gradient_test_1.grad_params)
    reference = np.array(
        [chi2_gradient_test_1.gradient[key] for key in
         chi2_gradient_test_1.grad_params])
    np.testing.assert_almost_equal(reference / result, 1., decimal=4)

    # changing something w/o updating: gives old values
    event.model.parameters.t_0 += 0.1
    result_1 = event.chi2_gradient
    np.testing.assert_almost_equal(result_1 / result, 1.)

    # run fit_fluxes: gives new values
    event.fit_fluxes()
    result_2 = event.calculate_chi2_gradient(chi2_gradient_test_1.grad_params)
    assert all(result != result_2)

    # Go back to previous results and test that calculate_chi2_gradient gives
    # results that don't match anything.
    event.model.parameters.t_0 -= 0.1
    result_3 = event.calculate_chi2_gradient(chi2_gradient_test_1.grad_params)
    assert all(result != result_3)
    assert all(result_2 != result_3)
    result_4 = event.get_chi2_gradient(chi2_gradient_test_1.grad_params)
    np.testing.assert_almost_equal(result_4 / result, 1.)


def test_chi2_gradient_2():
    """
    Double-up on the gradient test, i.e. duplicate the dataset and check that
    the gradient values double.
    """
    reference = np.array(
        [chi2_gradient_test_1.gradient[key] for key in
         chi2_gradient_test_1.grad_params])

    data = mm.MulensData(file_name=SAMPLE_FILE_02)
    event = mm.Event(
        datasets=[data, data],
        model=mm.Model(chi2_gradient_test_1.parameters),
        fix_blend_flux={data: 0.})
    result_0 = event.get_chi2_gradient(chi2_gradient_test_1.grad_params)
    result_1 = event.fits[1].chi2_gradient
    np.testing.assert_almost_equal(2. * reference / result_0, 1., decimal=4)
    np.testing.assert_almost_equal(2. * result_1 / result_0, 1.)
    np.testing.assert_almost_equal(reference / result_1, 1., decimal=4)

    # Change something and test fit.get_chi2_gradient
    event.model.parameters.t_0 += 0.1
    result_2 = event.fits[1].get_chi2_gradient(
        chi2_gradient_test_1.grad_params)
    assert result_2[0] != result_1[0]
    event.model.parameters.t_0 -= 0.1
    result_3 = event.fits[1].get_chi2_gradient(
        chi2_gradient_test_1.grad_params)
    np.testing.assert_almost_equal(result_3 / reference, 1., decimal=4)


def _test_event_chi2_gradient_rho():
    """
    test calculation of chi2 gradient including finite source effects
    MB08310 is used as an example
    """
    kwargs = {'comments': ["\\", "|"]}
    datasets = [
        mm.MulensData(file_name=SAMPLE_FILE_310_01, bandpass='R', **kwargs),
        mm.MulensData(file_name=SAMPLE_FILE_310_02, bandpass='U', **kwargs),
        mm.MulensData(file_name=SAMPLE_FILE_310_03, bandpass='I', **kwargs)]

    (gamma_I, gamma_V) = (0.44, 0.72)
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
    params = parameters.keys()
    gradient = {'t_0': 1283513.3068849628, 'u_0': 20492801.742886964,
                't_E': -9573.3589902395161, 'rho': -1503911.2409404013}
    reference = np.array([gradient[key] for key in params])

    # Create event and run test
    event = mm.Event(model=model, datasets=datasets)
    # result = event.get_chi2_gradient(list(params), fit_blending=False)

    # print(result)
    # print(reference)
    # np.testing.assert_almost_equal(reference / result, 1., decimal=2)


# ----------
# Event.get_ref_fluxes() Tests
def test_get_ref_fluxes():
    """Test Event.get_ref_fluxes()"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    event = mm.Event(data, model)

    # New (v2.0) method for fixing the blending
    (f_s_1, f_b_1) = event.get_ref_fluxes()
    event.fix_blend_flux[data] = 0.
    event.fit_fluxes()
    (f_s_2, f_b_2) = event.get_ref_fluxes()
    event.fix_blend_flux = {}
    event.fit_fluxes()
    (f_s_3, f_b_3) = event.get_ref_fluxes()

    assert f_b_2 == 0.
    assert f_s_1 == f_s_3
    assert f_b_1 == f_b_3
    np.testing.assert_almost_equal((f_s_1 + f_b_1)/f_s_2, 1., decimal=3)
    # Table 1 of Poleski et al. 2014:
    np.testing.assert_almost_equal(f_b_1 / f_s_1, 0.016, decimal=3)


def test_get_ref_fluxes_binary_source():
    """
    test for case with multiple datasets
    """
    # Setup
    (model, model_1, model_2) = generate_binary_source_models()
    (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)
    event = mm.Event([data_1, data_2], model)

    # default: ref = dataset 0
    (source_fluxes_1, blend_flux_1) = event.get_ref_fluxes()
    np.testing.assert_almost_equal(source_fluxes_1, np.array([100., 300.]))
    np.testing.assert_almost_equal(blend_flux_1, 50.)

    # set data_ref != dataset 0
    event.data_ref = 1
    (source_fluxes_2, blend_flux_2) = event.get_ref_fluxes()
    np.testing.assert_almost_equal(source_fluxes_2, np.array([20., 30.]))
    np.testing.assert_almost_equal(blend_flux_2, 50.)

    # set data_ref using a dataset
    event.data_ref = data_1
    (source_fluxes_1, blend_flux_1) = event.get_ref_fluxes()
    np.testing.assert_almost_equal(source_fluxes_1, np.array([100., 300.]))
    np.testing.assert_almost_equal(blend_flux_1, 50.)


class TestDataRef(unittest.TestCase):
    def test_1(self):
        """
        Try get_ref_fluxes() with an event with duplicated data
        """
        (model, model_1, model_2) = generate_binary_source_models()
        (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)
        event_1 = mm.Event([data_1, data_2, data_2], model)
        with self.assertRaises(ValueError):
            event_1.data_ref = data_2


# --------
# Event.get_flux_for_dataset() Tests
def test_get_flux_for_dataset():
    """
    Test that get_flux_for_dataset behaves as expected under
    various conditions.
    """

    # Setup
    (model, model_1, model_2) = generate_binary_source_models()
    (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)
    event = mm.Event([data_1, data_2], model)

    # first run: works
    (source_fluxes_init, blend_flux_init) = event.get_flux_for_dataset(1)
    np.testing.assert_almost_equal(source_fluxes_init, np.array([20., 30.]))
    np.testing.assert_almost_equal(blend_flux_init, 50.)
    (source_fluxes, blend_flux) = event.get_flux_for_dataset(data_2)
    np.testing.assert_almost_equal(source_fluxes, np.array([20., 30.]))
    np.testing.assert_almost_equal(blend_flux, 50.)

    # change model parameters w/o updating: gives old values
    event.model.parameters.t_E = 30.
    (source_fluxes, blend_flux) = event.get_flux_for_dataset(1)
    np.testing.assert_almost_equal(source_fluxes, np.array([20., 30.]))
    np.testing.assert_almost_equal(blend_flux, 50.)

    # run fit_fluxes(): gives new values
    event.fit_fluxes()
    (source_fluxes, blend_flux) = event.get_flux_for_dataset(1)
    assert source_fluxes[0] < source_fluxes_init[0]
    assert source_fluxes[1] < source_fluxes_init[1]
    assert blend_flux != blend_flux_init

    # change fit_blend w/o updating: gives old values
    event.fix_blend_flux = {data_2: 0.}
    (source_fluxes_fix, blend_flux_fix) = event.get_flux_for_dataset(1)
    np.testing.assert_almost_equal(blend_flux_fix, blend_flux)

    # run fit_fluxes(): gives new values
    event.fit_fluxes()
    (source_fluxes_fix, blend_flux_fix) = event.get_flux_for_dataset(1)
    np.testing.assert_almost_equal(blend_flux_fix, 0.)


def test_get_chi2_per_point():
    """
    test format of output: access a specific point in an event with multiple
    datasets
    """

    (model, model_1, model_2) = generate_binary_source_models()
    (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)

    # Modify data_2 to make two outliers:
    n = 100
    data_2.flux[n] += data_2.err_flux[n]
    data_2.flux[n + 1] -= data_2.err_flux[n + 1]
    chi2_exp_2 = np.zeros(len(data_2.time))
    chi2_exp_2[n:n + 2] = 1.

    # Create event and perform test
    event = mm.Event([data_1, data_2], model)
    np.testing.assert_almost_equal(
        event.get_chi2_per_point()[1], chi2_exp_2, decimal=6)
    np.testing.assert_almost_equal(
        event.get_chi2_per_point(bad=True)[1], chi2_exp_2, decimal=6)

    # Make tests for when bad data are present:
    bad_1 = data_1.bad.copy()
    bad_1[n] = True
    data_1.bad = bad_1
    bad_2 = data_2.bad.copy()
    bad_2[n] = True
    data_2.bad = bad_2
    chi2_exp_2_sumA = 0.99799855
    chi2_exp_2_sumB = 2.00200630

    chi2_per_point = event.get_chi2_per_point()
    assert np.isnan(np.sum(chi2_per_point[0]))
    np.testing.assert_almost_equal(np.nansum(chi2_per_point[0]), 0.)
    assert np.isnan(np.sum(chi2_per_point[1]))
    np.testing.assert_almost_equal(
        np.nansum(chi2_per_point[1]), chi2_exp_2_sumA)

    chi2_per_point = event.get_chi2_per_point(bad=True)
    np.testing.assert_almost_equal(np.sum(chi2_per_point[0]), 0.)
    np.testing.assert_almost_equal(np.sum(chi2_per_point[1]), chi2_exp_2_sumB)


# ---------
class TestFixedFluxes(unittest.TestCase):
    """test various combinations of fixed source and blend fluxes)"""

    def setUp(self):
        self.model = mm.Model({'t_0': 8000., 'u_0': 0.3, 't_E': 25.})
        self.f_source_1 = 1.
        self.f_blend_1 = 3.
        self.f_source_2 = 4.
        self.f_blend_2 = 1.5
        self.dt = 2.0
        self.n_tE = 100
        self.flux_uncertainty = 0.01

        self.constrain_1 = 1.
        self.constrain_2 = 3.
        self.constrain_3 = 1.1

        self.expected_1 = 5.455106739005208
        self.expected_2 = 4.485035579668404
        self.expected_3 = 3.910213477978024
        self.expected_4 = 2.898667067448723
        self.expected_5 = 2.544893260994792
        self.expected_6 = 2.5133293255121805

        self.generate_two_fake_datasets()

    def generate_two_fake_datasets(self):
        """
        Generate two perfect datasets with different source and blend fluxes
        """
        t_0 = self.model.parameters.t_0
        delta_t = self.n_tE * self.model.parameters.t_E
        times = np.arange(t_0-delta_t, t_0+delta_t, self.dt)

        self.data_1 = self.generate_dataset(
            self.f_source_1, self.f_blend_1, times)

        self.data_2 = self.generate_dataset(
            self.f_source_2, self.f_blend_2, times+self.dt/2.)

        self.datasets = [self.data_1, self.data_2]

    def generate_dataset(self, f_source, f_blend, times):
        """generate a single dataset without noise"""
        flux = f_source * self.model.get_magnification(times) + f_blend
        err = np.zeros(len(times)) + self.flux_uncertainty
        data = mm.MulensData([times, flux, err], phot_fmt='flux')
        return data

    def extract_fluxes(self, event):
        event.fit_fluxes()
        fluxes_1 = event.get_flux_for_dataset(self.data_1)
        fluxes_2 = event.get_flux_for_dataset(self.data_2)
        return (fluxes_1, fluxes_2)

    def test_fixed_blend_fluxes_1(self):
        """
        test fix no blending flux for second dataset and all other fluxes free
        """
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_2: 0.})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(fluxes_1[0][0], self.f_source_1)
        np.testing.assert_almost_equal(fluxes_1[1], self.f_blend_1)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.expected_1)
        np.testing.assert_almost_equal(fluxes_2[1], 0.)

    def test_fixed_blend_fluxes_2(self):
        """test one is fixed at value, the other is free"""
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_2: self.constrain_1})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(fluxes_1[0][0], self.f_source_1)
        np.testing.assert_almost_equal(fluxes_1[1], self.f_blend_1)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.expected_2)
        np.testing.assert_almost_equal(fluxes_2[1], 1.)

    def test_fixed_blend_fluxes_3(self):
        """test both fixed at zero"""
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_1: 0., self.data_2: 0.})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(fluxes_1[0][0], self.expected_3)
        np.testing.assert_almost_equal(fluxes_1[1], 0.)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.expected_1)
        np.testing.assert_almost_equal(fluxes_2[1], 0.)

    def test_fixed_source_and_blend_fluxes(self):
        """test one blend fixed at zero and the other source fixed"""
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={self.data_1: self.constrain_3},
            fix_blend_flux={self.data_2: self.constrain_2})
        event.fit_fluxes()
        (fluxes_1, fluxes_2) = event.fluxes

        np.testing.assert_almost_equal(fluxes_1[0], self.constrain_3)
        np.testing.assert_almost_equal(fluxes_1[1], self.expected_4)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.expected_5)
        np.testing.assert_almost_equal(fluxes_2[1], self.constrain_2)

    def test_fixed_source_fluxes_1(self):
        """test one is fixed at value, the other is free"""
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={self.data_2: self.constrain_2})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            (fluxes_1[0][0], fluxes_1[1]), (self.f_source_1, self.f_blend_1))
        np.testing.assert_almost_equal(fluxes_2[0][0], self.constrain_2)

        np.testing.assert_almost_equal(fluxes_2[1], self.expected_6)

    def test_fixed_source_fluxes_2(self):
        """test both fixed at values"""
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={
                self.data_2: self.constrain_2,
                self.data_1: self.constrain_3})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(fluxes_1[0][0], self.constrain_3)
        np.testing.assert_almost_equal(fluxes_1[1], self.expected_4)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.constrain_2)
        np.testing.assert_almost_equal(fluxes_2[1], self.expected_6)

    def test_free_fluxes(self):
        """test fit for all fluxes free"""
        event = mm.Event(datasets=self.datasets, model=self.model)
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(fluxes_1[0][0], self.f_source_1)
        np.testing.assert_almost_equal(fluxes_1[1], self.f_blend_1)
        np.testing.assert_almost_equal(fluxes_2[0][0], self.f_source_2)
        np.testing.assert_almost_equal(fluxes_2[1], self.f_blend_2)


# ---------
class TestFixedFluxRatios(unittest.TestCase):
    """test various combinations of fixed flux ratios"""

    def setUp(self):
        # Flux ratios
        self.q_I = 0.01
        self.q_V = 0.002

        self.gamma = {'I': 0.44, 'V': 0.72}

        # Model based on binary source
        t_0_2 = 8001.
        self.model = mm.Model(
            {'t_0_1': 8000., 'u_0_1': 0.3, 't_0_2': t_0_2, 'u_0_2': 0.001,
             't_E': 25., 'rho_2': 0.002})
        d_t = 3. * self.model.parameters.t_star_2
        self.model.set_magnification_methods(
            [t_0_2 - d_t, 'finite_source_LD_Yoo04', t_0_2 + d_t], source=2)
        self.model.set_limb_coeff_gamma('I', self.gamma['I'])
        self.model.set_limb_coeff_gamma('V', self.gamma['V'])

        self.generate_fake_datasets()

    def generate_fake_datasets(self):
        """
        Generate perfect datasets with different source and blend fluxes
        """
        model_1 = mm.Model(self.model.parameters.source_1_parameters)
        model_2 = mm.Model(self.model.parameters.source_2_parameters)
        model_2.set_magnification_methods(self.model._methods[2])
        model_2.set_limb_coeff_gamma('I', self.gamma['I'])
        model_2.set_limb_coeff_gamma('V', self.gamma['V'])

        def gen_data(
                f_source_1, f_blend, q_flux, times, bandpass=None, **kwargs):
            """generate perfect data for a given set of properties"""

            if bandpass == 'I':
                gamma = self.gamma['I']
            elif bandpass == 'V':
                gamma = self.gamma['V']
            else:
                gamma = None

            mag_1 = model_1.get_magnification(times)
            mag_2 = model_2.get_magnification(times, gamma=gamma)
            f_source_2 = f_source_1 * q_flux
            flux = f_source_1 * mag_1 + f_source_2 * mag_2 + f_blend
            err = np.zeros(len(times)) + 0.01
            data = mm.MulensData(
                [times, flux, err], phot_fmt='flux', bandpass=bandpass,
                **kwargs)

            return data

        def add_data(properties, label=None):
            """create data in two bands for each fake observatory"""
            times = np.arange(
                properties['t_start'], properties['t_stop'],
                properties['dt_I'])
            data_I = gen_data(
                properties['f_source_I'], properties['f_blend_I'], self.q_I,
                times, bandpass='I',
                plot_properties={'label': '{0} I'.format(label)})
            self.datasets.append(data_I)
            times = np.arange(
                properties['t_start'], properties['t_stop'],
                properties['dt_V'])
            data_V = gen_data(
                properties['f_source_V'], properties['f_blend_V'], self.q_V,
                times, bandpass='V',
                plot_properties={'label': '{0} V'.format(label)})
            self.datasets.append(data_V)

        self.datasets = []
        self.expected_fluxes = []
        n_tE = 10
        self.data_properties = {
            # dense "Survey" data
            'survey_1': {
                't_start': (self.model.parameters.t_0_1 -
                            n_tE * self.model.parameters.t_E),
                't_stop': (self.model.parameters.t_0_1 +
                           n_tE * self.model.parameters.t_E),
                'dt_I': 0.04, 'dt_V': 0.4, 'f_source_I': 1., 'f_blend_I': 0.2,
                'f_source_V': 0.8, 'f_blend_V': 0.3},
            # sparse "Survey" data
            'survey_2': {
                't_start': (self.model.parameters.t_0_1 -
                            n_tE * self.model.parameters.t_E+1.01),
                't_stop': (self.model.parameters.t_0_1 +
                           n_tE * self.model.parameters.t_E+1.01),
                'dt_I': 1.0, 'dt_V': 2.0, 'f_source_I': 1.1, 'f_blend_I': 0.15,
                'f_source_V': 0.81, 'f_blend_V': 0.2},
            # "FollowUp" data
            'followup_1': {
                't_start': (self.model.parameters.t_0_1 - 1.1),
                't_stop': (self.model.parameters.t_0_2 +
                           self.model.parameters.t_E),
                'dt_I': 0.001, 'dt_V': 0.01, 'f_source_I': 2.3,
                'f_blend_I': 0.5, 'f_source_V': 1.8, 'f_blend_V': 0.65}}

        # The data_keys list is intentional to ensure that survey_1, I is the
        # first (and reference) dataset
        self.data_keys = ['survey_1', 'survey_2', 'followup_1']
        for key in self.data_keys:
            add_data(self.data_properties[key], label=key)
            self.expected_fluxes.append(
                ([self.data_properties[key]['f_source_I'],
                  self.q_I * self.data_properties[key]['f_source_I']],
                 self.data_properties[key]['f_blend_I']))
            self.expected_fluxes.append(
                ([self.data_properties[key]['f_source_V'],
                  self.q_V * self.data_properties[key]['f_source_V']],
                 self.data_properties[key]['f_blend_V']))

    def extract_fluxes(self, event):
        """ Return a list of fitted fluxes for all datasets"""
        event.fit_fluxes()
        fluxes = []
        for (i, dataset) in enumerate(self.datasets):
            fluxes.append(event.get_flux_for_dataset(i))

        return fluxes

    def test_both_q_fixed(self):
        """test the case that q_I is fixed"""
        q_values = {'I': 0.013, 'V': 0.005}

        event = mm.Event(datasets=self.datasets, model=self.model)
        event.fix_source_flux_ratio = q_values
        fluxes = self.extract_fluxes(event)
        for (i, dataset) in enumerate(self.datasets):
            np.testing.assert_almost_equal(
                fluxes[i][0][1] / fluxes[i][0][0], q_values[dataset.bandpass])

            assert event.get_chi2_for_dataset(i) > 1

    def test_fixed_q_I_flux(self):
        """test the case that q_I is fixed"""
        q_I_value = 0.012

        event = mm.Event(datasets=self.datasets, model=self.model)
        event.fix_source_flux_ratio = {'I': q_I_value}
        fluxes = self.extract_fluxes(event)
        for (i, dataset) in enumerate(self.datasets):
            if dataset.bandpass == 'I':
                # the ratio of q_I should be identical to the set value
                np.testing.assert_almost_equal(
                    fluxes[i][0][1] / fluxes[i][0][0], q_I_value)
                assert event.get_chi2_for_dataset(i) > 1
            elif dataset.bandpass == 'V':
                # the ratio of q_V should be the input value
                np.testing.assert_almost_equal(
                    fluxes[i][0][1] / fluxes[i][0][0], self.q_V)
                np.testing.assert_almost_equal(
                    event.get_chi2_for_dataset(i), 0.)

    def test_free_fluxes(self):
        """test both q_flux free. Should give original values"""
        event = mm.Event(datasets=self.datasets, model=self.model)
        fluxes = self.extract_fluxes(event)
        for (i, dataset) in enumerate(self.datasets):
            np.testing.assert_almost_equal(
                fluxes[i][0][0] / self.expected_fluxes[i][0][0], 1.)
            np.testing.assert_almost_equal(
                fluxes[i][0][1] / self.expected_fluxes[i][0][1], 1.)
            np.testing.assert_almost_equal(
                fluxes[i][1] / self.expected_fluxes[i][1], 1.)
            np.testing.assert_almost_equal(event.get_chi2_for_dataset(i), 0.)

# Tests to add:
#
# properties: coords, model, datasets, data_ref, sum_function?
