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
    # assert ev.fit is None or isinstance(ev.fit, mm.Fit)

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 427.20382, decimal=4,
                                   err_msg='problem in resulting chi2')

    ev.sum_function = 'numpy.sum'  # JCY what is the purpose of this test?
    np.testing.assert_almost_equal(
        ev.get_chi2(), 427.20382, decimal=4, err_msg='problem with numpy.sum')

    # Old method of fixing the blending
    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')

    # New method of fixing the blending
    ev.fix_blend_flux = {}
    np.testing.assert_almost_equal(ev.get_chi2(), 427.20382, decimal=4)

    fix_blend_flux = {}
    for dataset in ev.datasets:
        fix_blend_flux[dataset] = 0.

    ev.fix_blend_flux = fix_blend_flux
    np.testing.assert_almost_equal(ev.get_chi2(), 459.09826, decimal=4)


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
    mod.set_datasets([data_1, data_2])
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

    # Old method of fixing the blending
    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 2.*459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')

    # New method of fixing the blending
    # Both datasets have zero blending
    ev.fix_blend_flux = {data_1: 0., data_2: 0.}
    chi2_no_blend_2 = ev.get_chi2()
    assert isinstance(chi2_no_blend_2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(
        float(chi2_no_blend), 2.*459.09826, decimal=4,
        err_msg='problem in resulting chi2 for fixed no blending')


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

    # Old method for fixing blending
    # for test in [test_1]:  # , test_2]:
    #     (parameters, params, gradient) = test
    #     event = mm.Event(model=mm.Model(parameters), **kwargs)
    #     event.fit_fluxes()
    #     result = event.chi2_gradient(params, fit_blending=False)
    #
    #     reference = np.array([gradient[key] for key in params])
    #     np.testing.assert_almost_equal(reference/result, 1., decimal=1)

    # New method for fixing blending
    for test in [chi2_gradient_test_1]:  # , chi2_gradient_test_2]:
        event = mm.Event(
            model=mm.Model(test.parameters), fix_blend_flux={data: 0.}, **kwargs)
        result = event.get_chi2_gradient(test.grad_params)
        reference = np.array([test.gradient[key] for key in test.grad_params])
        np.testing.assert_almost_equal(reference/result, 1., decimal=4)

        result = event.chi2_gradient
        np.testing.assert_almost_equal(reference/result, 1., decimal=4)

class TestGradient(unittest.TestCase):

    def test_gradient_init_1(self):
        """test that fit_fluxes/update must be run in order to update the
        gradient"""
        # before fit_fluxes --> error
        with self.assertRaises(TypeError):
            data = mm.MulensData(file_name=SAMPLE_FILE_02)
            event = mm.Event(
                datasets=[data],
                model=mm.Model(chi2_gradient_test_1.parameters),
                fix_blend_flux={data: 0.})
            event.chi2_gradient()

def test_chi2_gradient():
    """test that fit_fluxes/update must be run in order to update the
    gradient"""
    data = mm.MulensData(file_name=SAMPLE_FILE_02)
    event = mm.Event(
        datasets=[data],
        model=mm.Model(chi2_gradient_test_1.parameters),
        fix_blend_flux={data: 0.})
    result = event.get_chi2_gradient(chi2_gradient_test_1.grad_params)
    reference = np.array(
        [chi2_gradient_test_1.gradient[key] for key in
         chi2_gradient_test_1.grad_params])
    np.testing.assert_almost_equal(reference / result, 1., decimal=1)

    # changing something w/o updating: gives old values
    event.model.parameters.t_0 += 0.1
    result_1 = event.chi2_gradient
    np.testing.assert_almost_equal(result_1 / result, 1.)

    def is_different(result, new_result):
        assert result[0] != new_result[0]
        assert result[1] != new_result[1]
        assert result[2] != new_result[2]

    # run fit_fluxes: gives new values
    event.fit_fluxes()
    result_2 = event.calc_chi2_gradient(chi2_gradient_test_1.grad_params)
    is_different(result, result_2)

    # Go back to previous results and test that calc_chi2_gradient gives
    # results that don't match anything.
    event.model.parameters.t_0 -= 0.1
    result_3 = event.calc_chi2_gradient(chi2_gradient_test_1.grad_params)
    is_different(result, result_3)
    is_different(result_2, result_3)
    result_4 = event.get_chi2_gradient(chi2_gradient_test_1.grad_params)
    np.testing.assert_almost_equal(result_4 / result, 1.)

def test_chi2_gradient_2():
    # Double-up on the gradient test, i.e. duplicate the dataset and check that the
    # gradient values double.
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
    result_2 = event.fits[1].get_chi2_gradient(chi2_gradient_test_1.grad_params)
    assert result_2[0] != result_1[0]
    event.model.parameters.t_0 -= 0.1
    result_3 = event.fits[1].get_chi2_gradient(chi2_gradient_test_1.grad_params)
    np.testing.assert_almost_equal(result_3 / reference, 1., decimal=4)

def test_get_ref_fluxes():
    """Test Event.get_ref_fluxes()"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    event = mm.Event(data, model)

    # Old method for fixing the blending
    (f_s_1, f_b_1) = event.get_ref_fluxes()
    (f_s_2, f_b_2) = event.get_ref_fluxes(fit_blending=False)
    (f_s_3, f_b_3) = event.get_ref_fluxes(fit_blending=True)

    assert f_b_2 == 0.
    assert f_s_1 == f_s_3
    assert f_b_1 == f_b_3
    np.testing.assert_almost_equal((f_s_1 + f_b_1)/f_s_2, 1., decimal=3)
    # Table 1 of Poleski et al. 2014:
    np.testing.assert_almost_equal(f_b_1 / f_s_1, 0.016, decimal=3)

    # New method for fixing the blending
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

    # Old method of setting flux ratios
    # model.set_source_flux_ratio(1.)
    # np.testing.assert_almost_equal(model.magnification(time), (mag_1+mag_2)/2.)
    # assert event.get_chi2() > 1.
    # model.set_source_flux_ratio(3.)
    # np.testing.assert_almost_equal(event.get_chi2(), 0.)

    # New method of setting flux ratios
    event.fix_q_flux = {data: 1.}
    assert event.get_chi2() > 1.
    event.fix_q_flux = {data: 3.}
    np.testing.assert_almost_equal(event.get_chi2(), 0.)


def generate_binary_source_models():
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.05,
        't_0_2': 5100., 'u_0_2': 0.15, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    return (model, model_1, model_2)

def generate_binary_source_datasets(model_1, model_2):
    # prepare fake data:
    time = np.linspace(4900., 5200., 600)
    mag_1 = model_1.magnification(time)
    mag_2 = model_2.magnification(time)
    flux = 100. * mag_1 + 300. * mag_2 + 50.
    data_1 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')
    flux = 20. * mag_1 + 30. * mag_2 + 50.
    data_2 = mm.MulensData(data_list=[time, flux, 1.+0.*time], phot_fmt='flux')

    return (data_1, data_2)

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

    # Old method of fixing flux ratios
    # data_1.bandpass = 'some'
    # event.model.set_source_flux_ratio_for_band('some', 3.)
    # np.testing.assert_almost_equal(event.get_chi2_for_dataset(0), 0.)

    # New method of fixing flux ratios
    data_1.bandpass = 'some'
    event.fix_q_flux['some'] = 3.
    event.fit_fluxes()
    np.testing.assert_almost_equal(event.get_chi2_for_dataset(0), 0.)

def test_get_flux_for_dataset():
    """ Test that get_flux_for_dataset behaves as expected under various
    conditions. """

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


def test_get_ref_fluxes():
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
    event.data_ref = 1 # Does not work
    (source_fluxes_2, blend_flux_2) = event.get_ref_fluxes()
    np.testing.assert_almost_equal(source_fluxes_2, np.array([20., 30.]))
    np.testing.assert_almost_equal(blend_flux_2, 50.)

    # set data_ref using a dataset
    event.data_ref = data_1 # Does not work
    (source_fluxes_1, blend_flux_1) = event.get_ref_fluxes()
    np.testing.assert_almost_equal(source_fluxes_1, np.array([100., 300.]))
    np.testing.assert_almost_equal(blend_flux_1, 50.)

# Try get_ref_fluxes() with an event with duplicated data:
class TestDataRef(unittest.TestCase):
    def test_1(self):
        (model, model_1, model_2) = generate_binary_source_models()
        (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)
        with self.assertRaises(ValueError):
            event_1 = mm.Event([data_1, data_2, data_2], model)
            event_1.data_ref = data_2

def test_get_chi2_per_point():
    """test format of output: access a specific point in an event with multiple
    datasets"""

    # Setup
    (model, model_1, model_2) = generate_binary_source_models()
    (data_1, data_2) = generate_binary_source_datasets(model_1, model_2)

    # Modify data_2 to make two outliers:
    n = 100
    data_2.flux[n] += data_2.err_flux[n]
    data_2.flux[n + 1] -= data_2.err_flux[n + 1]
    chi2_exp = np.zeros(len(data_2.time))
    chi2_exp[n:n + 2] = 1.

    # Create event and perform test
    event = mm.Event([data_1, data_2], model)
    np.testing.assert_almost_equal(
        event.get_chi2_per_point()[1], chi2_exp, decimal=6)


class TestFixedFluxes(unittest.TestCase):
    """test various combinations of fixed source and blend fluxes)"""

    def setUp(self):
        self.model = mm.Model(
            {'t_0': 8000., 'u_0': 0.3, 't_E': 25.})
        self.generate_two_fake_datasets()

    def generate_two_fake_datasets(self):
        """
        Generate two perfect datasets with different source and blend fluxes
        """
        self.f_source_1, self.f_blend_1 = 1.0, 3.0
        self.f_source_2, self.f_blend_2 = 4.0, 1.5

        def gen_data(f_source, f_blend, times):
            flux = f_source * self.model.magnification(times) + f_blend
            err = np.zeros(len(times)) + 0.01
            data = mm.MulensData([times, flux, err], phot_fmt='flux')
            return data

        dt = 2.0
        n = 100
        self.data_1 = gen_data(
            self.f_source_1, self.f_blend_1,
            np.arange(self.model.parameters.t_0 - n * self.model.parameters.t_E,
                      self.model.parameters.t_0 + n * self.model.parameters.t_E,
                      dt))
        self.data_2 = gen_data(
            self.f_source_2, self.f_blend_2,
            np.arange(
                self.model.parameters.t_0 - n * self.model.parameters.t_E + dt / 2.,
                self.model.parameters.t_0 + n * self.model.parameters.t_E + dt / 2.,
                      dt))

        self.datasets = [self.data_1, self.data_2]

    def extract_fluxes(self, event):
        event.fit_fluxes()
        fluxes_1 = event.get_flux_for_dataset(self.data_1)
        fluxes_2 = event.get_flux_for_dataset(self.data_2)
        print(fluxes_1)
        print(fluxes_2)
        return (fluxes_1, fluxes_2)

    def test_free_fluxes(self):
        # test both free
        event = mm.Event(datasets=self.datasets, model=self.model)
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1, (self.f_source_1, self.f_blend_1))
        np.testing.assert_almost_equal(
            fluxes_2, (self.f_source_2, self.f_blend_2))

    def test_fixed_blend_fluxes_1(self):
        # test one is fixed at zero, the other is free.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_2: 0.})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1, (self.f_source_1, self.f_blend_1))
        np.testing.assert_almost_equal(
            fluxes_2[0] / (self.f_source_2 + self.f_blend_2), 1., decimal=2)
        np.testing.assert_almost_equal(
            fluxes_2[1], 0.)

    def test_fixed_blend_fluxes_2(self):
        # test one is fixed at value, the other is free.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_2: 1.})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1, (self.f_source_1, self.f_blend_1))
        np.testing.assert_almost_equal(
            fluxes_2[0] / (self.f_source_2 + self.f_blend_2 - 1.), 1., decimal=2)
        np.testing.assert_almost_equal(
            fluxes_2[1], 1.)

    def test_fixed_blend_fluxes_3(self):
        # test both fixed at zero.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_blend_flux={self.data_1: 0., self.data_2: 0.})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1[0] / (self.f_source_1 + self.f_blend_1), 1., decimal=1)
        np.testing.assert_almost_equal(
            fluxes_1[1], 0.)
        np.testing.assert_almost_equal(
            fluxes_2[0] / (self.f_source_2 + self.f_blend_2), 1., decimal=2)
        np.testing.assert_almost_equal(
            fluxes_2[1], 0.)

    def test_fixed_source_fluxes_1(self):
        # test one is fixed at value, the other is free.
        fixed_source_flux = 3.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={self.data_2: fixed_source_flux})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1, (self.f_source_1, self.f_blend_1))
        np.testing.assert_almost_equal(
            fluxes_2[0], fixed_source_flux)
        np.testing.assert_almost_equal(
            fluxes_2[1] / (self.f_source_2 - fixed_source_flux + self.f_blend_2),
            1., decimal=2)

    def test_fixed_source_fluxes_2(self):
        # test both fixed at values.
        fixed_source_flux_1 = 1.1
        fixed_source_flux_2 = 3.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={
                self.data_2: fixed_source_flux_2,
                self.data_1: fixed_source_flux_1})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1[0], fixed_source_flux_1)
        np.testing.assert_almost_equal(
            fluxes_1[1] / (self.f_source_1 - fixed_source_flux_1 + self.f_blend_1),
            1., decimal=2)
        np.testing.assert_almost_equal(
            fluxes_2[0], fixed_source_flux_2)
        np.testing.assert_almost_equal(
            fluxes_2[1] / (self.f_source_2 - fixed_source_flux_2 + self.f_blend_2),
            1., decimal=2)


    def test_fixed_source_and_blend_fluxes(self):
        # test one blend fixed at zero and the other source fixed.
        fixed_source_flux_1 = 1.1
        fixed_blend_flux_2 = 3.
        event = mm.Event(
            datasets=self.datasets, model=self.model,
            fix_source_flux={self.data_1: fixed_source_flux_1},
            fix_blend_flux={self.data_2: fixed_blend_flux_2})
        (fluxes_1, fluxes_2) = self.extract_fluxes(event)

        np.testing.assert_almost_equal(
            fluxes_1[0], fixed_source_flux_1)
        np.testing.assert_almost_equal(
            fluxes_1[1] / (self.f_source_1 - fixed_source_flux_1 + self.f_blend_1),
            1., decimal=2)
        np.testing.assert_almost_equal(
            fluxes_2[0] / (self.f_source_2 + self.f_blend_2 - fixed_blend_flux_2),
            1, decimal=1)
        np.testing.assert_almost_equal(
            fluxes_2[1], fixed_blend_flux_2)

# Tests to add:
#
#     fix_q_flux: in the case that q_I is fixed, e.g. for two datasets with
#       band=I and one with band=V

# test chi2 vs get_chi2:
#      get_chi2 always updates after something changes in the model, chi2 does
#        not.
#
# properties: coords, model, datasets, data_ref, sum_function?