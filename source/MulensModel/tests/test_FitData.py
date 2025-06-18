import ast

import numpy as np
from numpy.testing import assert_almost_equal as almost
from numpy.testing import assert_allclose
import pytest
import unittest
from os.path import join

import MulensModel as mm
from fortran_files import FortranSFitFile51

dir_1 = join(mm.DATA_PATH, 'photometry_files', 'OB140939')
dir_2 = join(mm.DATA_PATH, 'unit_test_files')
dir_3 = join(mm.DATA_PATH, 'ephemeris_files')
dir_4 = join(dir_2, 'fspl_derivs')

SAMPLE_FILE_02 = join(dir_1, 'ob140939_OGLE.dat')  # HJD'
SAMPLE_FILE_02_REF = join(dir_2, 'ob140939_OGLE_ref_v2.dat')  # HJD'
SAMPLE_FILE_03 = join(dir_1, 'ob140939_Spitzer.dat')  # HJD'
SAMPLE_FILE_03_EPH = join(dir_3, 'Spitzer_ephemeris_01.dat')  # UTC
SAMPLE_FILE_03_REF = join(dir_2, 'ob140939_Spitzer_ref_v2.dat')  # HJD'
SAMPLE_FILE_04_WF = join(mm.DATA_PATH, 'WFIRST_1827.dat')
SAMPLE_FILE_FSPL_51 = join(dir_4, 'fort.51')
SAMPLE_FILE_FSPL_61 = join(dir_4, 'fort.61')


# Note: default precision for assert_almost_equal (aka almost) is decimal = 7


def generate_model():
    """
    returns a model, time array, and magnification
    """

    # Create a PSPL model
    t_0 = 3583.
    u_0 = 0.3
    t_E = 12.

    t = np.linspace(t_0 - 3. * t_E, t_0 + 3. * t_E, 1000)
    pspl = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    A = pspl.get_magnification(t)

    return (pspl, t, A)


def generate_binary_model():
    """
    returns a binary source model, time array, and the magnification of
    both sources
    """

    # retrieve model 1
    (model_1, t, A_1) = generate_model()
    t_0_1 = model_1.parameters.t_0
    u_0_1 = model_1.parameters.u_0

    # create second model
    t_0_2 = 3570.
    u_0_2 = 0.25
    t_E = model_1.parameters.t_E

    model_2 = mm.Model({'t_0': t_0_2, 'u_0': u_0_2, 't_E': t_E})

    A_2 = model_2.get_magnification(t)

    # create separate binary model
    params = {'t_0_1': t_0_1, 'u_0_1': u_0_1, 't_0_2': t_0_2, 'u_0_2': u_0_2,
              't_E': t_E}
    binary_model = mm.Model(params)

    return (binary_model, t, A_1, A_2)


def generate_dataset(f_mod, t):
    """
    pass in f_mod and t, returns a MulensData
    """

    # error in measurement
    err = f_mod * 0.01

    my_dataset = mm.MulensData(data_list=[t, f_mod, err], phot_fmt='flux')

    return my_dataset


class TestSingleSourceFluxes(unittest.TestCase):

    def setUp(self):
        self.pspl, self.t, self.magnification = generate_model()

        # secrets
        self.f_s = 1.0
        self.f_b = 0.5
        # generate f_mod
        self.f_mod = self.f_s * self.magnification + self.f_b

        self.my_dataset = generate_dataset(self.f_mod, self.t)

    def _run_true_values_test(
            self, fix_source_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.pspl, dataset=self.my_dataset,
            fix_blend_flux=fix_blend_flux, fix_source_flux=fix_source_flux)
        my_fit.fit_fluxes()

        almost(my_fit.blend_flux, self.f_b)
        almost(my_fit.source_flux, self.f_s)

        # Test get_model_fluxes() for 1 source
        peak_index = 500
        mod_fluxes = my_fit.get_model_fluxes()
        almost(mod_fluxes[peak_index], self.my_dataset.flux[peak_index])

    def _run_arbitrary_values_test(
            self, fix_source_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.pspl, dataset=self.my_dataset,
            fix_blend_flux=fix_blend_flux, fix_source_flux=fix_source_flux)
        my_fit.fit_fluxes()

        if fix_blend_flux is not False:
            almost(my_fit.blend_flux, fix_blend_flux)

        if fix_source_flux is not False:
            almost(my_fit.source_flux, fix_source_flux)

    def test_all_free(self):
        self._run_true_values_test()
        self._run_true_values_test(fix_source_flux=False, fix_blend_flux=False)

    def test_fixed_blending(self):
        self._run_true_values_test(fix_blend_flux=0.5)
        self._run_arbitrary_values_test(fix_blend_flux=0.)
        self._run_arbitrary_values_test(fix_blend_flux=0.23)
        self._run_arbitrary_values_test(fix_blend_flux=-0.23)

    def test_fixed_source_flux(self):
        self._run_true_values_test(fix_source_flux=1.0)
        self._run_arbitrary_values_test(fix_source_flux=0.)
        self._run_arbitrary_values_test(fix_source_flux=0.23)
        self._run_arbitrary_values_test(fix_source_flux=-0.23)

    def test_all_fixed(self):
        self._run_true_values_test(fix_source_flux=1.0, fix_blend_flux=0.5)
        self._run_arbitrary_values_test(
            fix_source_flux=1.2, fix_blend_flux=-0.5)
        self._run_arbitrary_values_test(fix_source_flux=1.7, fix_blend_flux=0.)

    def test_data_magnification(self):
        dataset = self.my_dataset.copy()
        dataset.bad = (dataset.time < self.pspl.parameters.t_0)
        magnification = np.zeros(len(dataset.time))
        magnification[dataset.good] = self.pspl.get_magnification(
            dataset.time[dataset.good])

        my_fit = mm.FitData(model=self.pspl, dataset=dataset)
        assert my_fit.data_magnification is None
        my_fit.get_data_magnification()
        np.testing.assert_equal(
            my_fit.data_magnification, magnification)


class TestBinarySourceFluxes(unittest.TestCase):

    def setUp(self):
        self.f_s_1 = 1
        self.f_s_2 = 1.2
        self.f_b = 0.5

        self.model, self.t, self.A_1, self.A_2 = generate_binary_model()
        f_mod = self.f_s_1 * self.A_1 + self.f_s_2 * self.A_2 + self.f_b
        self.dataset = generate_dataset(f_mod, self.t)

    def _run_true_value_test(
            self, fix_source_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_source_flux=fix_source_flux, fix_blend_flux=fix_blend_flux)
        my_fit.update()
        almost(my_fit.blend_flux, self.f_b)
        almost(my_fit.source_fluxes[0], self.f_s_1)
        almost(my_fit.source_fluxes[1], self.f_s_2)
        assert isinstance(my_fit.source_fluxes, (np.ndarray))

        peak_index = 500
        mod_fluxes = my_fit.get_model_fluxes()
        almost(mod_fluxes[peak_index], self.dataset.flux[peak_index])

        assert (my_fit.chi2_per_point.shape == (self.dataset.n_epochs, ))

    def _run_arbitrary_value_test(
            self, fix_source_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_source_flux=fix_source_flux, fix_blend_flux=fix_blend_flux)
        my_fit.fit_fluxes()

        if fix_blend_flux is not False:
            almost(my_fit.blend_flux, fix_blend_flux)

        if fix_source_flux is not False:
            if fix_source_flux[0] is not False:
                almost(my_fit.source_fluxes[0], fix_source_flux[0])

            if fix_source_flux[1] is not False:
                almost(my_fit.source_fluxes[1], fix_source_flux[1])

    def _run_q_flux_test_true(self, fix_q_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_source_flux_ratio=fix_q_flux, fix_blend_flux=fix_blend_flux)
        my_fit.fit_fluxes()

        almost(my_fit.blend_flux, self.f_b)
        almost(my_fit.source_fluxes[0], self.f_s_1)
        almost(my_fit.source_fluxes[1], self.f_s_2)

    def _run_q_flux_test_arbitrary(
            self, fix_q_flux=False, fix_blend_flux=False):
        my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_source_flux_ratio=fix_q_flux, fix_blend_flux=fix_blend_flux)
        my_fit.fit_fluxes()

        if fix_blend_flux is not False:
            almost(my_fit.blend_flux, fix_blend_flux)

        almost(
            my_fit.source_fluxes[1]/my_fit.source_fluxes[0],
            fix_q_flux)
        assert isinstance(my_fit.source_fluxes, (np.ndarray))

    def test_value_error(self):
        with self.assertRaises(ValueError):
            self._run_true_value_test(fix_source_flux=1.0)

    def test_all_free(self):
        self._run_true_value_test()
        self._run_true_value_test(
            fix_source_flux=[False, False], fix_blend_flux=False)

    def test_fixed_source_true(self):
        self._run_true_value_test(
            fix_source_flux=[1., False], fix_blend_flux=False)
        self._run_true_value_test(
            fix_source_flux=[False, 1.2], fix_blend_flux=False)
        self._run_true_value_test(
            fix_source_flux=[1., 1.2], fix_blend_flux=False)

    def test_fixed_blend_true(self):
        self._run_true_value_test(fix_blend_flux=0.5)
        self._run_true_value_test(
            fix_source_flux=[1., False], fix_blend_flux=0.5)
        self._run_true_value_test(
            fix_source_flux=[False, 1.2], fix_blend_flux=0.5)

    def test_all_fixed_true(self):
        self._run_true_value_test(
            fix_source_flux=[1., 1.2], fix_blend_flux=0.5)

    def test_fixed_source_arbitrary(self):
        self._run_arbitrary_value_test(
            fix_source_flux=[1.2, False], fix_blend_flux=False)
        self._run_arbitrary_value_test(
            fix_source_flux=[False, 0.53], fix_blend_flux=False)
        self._run_arbitrary_value_test(
            fix_source_flux=[4.5, 0.67], fix_blend_flux=False)

    def test_fixed_blend_arbitrary(self):
        self._run_arbitrary_value_test(fix_blend_flux=0.)
        self._run_arbitrary_value_test(fix_blend_flux=2.3)
        self._run_arbitrary_value_test(fix_blend_flux=-0.5)
        self._run_arbitrary_value_test(
            fix_source_flux=[1.2, False], fix_blend_flux=0.78)
        self._run_arbitrary_value_test(
            fix_source_flux=[False, 0.53], fix_blend_flux=0.23)

    def test_all_fixed_arbitrary(self):
        self._run_arbitrary_value_test(
            fix_source_flux=[2.3, 0.45], fix_blend_flux=0.67)
        self._run_arbitrary_value_test(
            fix_source_flux=[2.3, 0.45], fix_blend_flux=0.)

    def test_q_flux_fixed(self):
        self._run_q_flux_test_true()
        self._run_q_flux_test_true(fix_q_flux=1.2)
        self._run_q_flux_test_true(fix_q_flux=1.2, fix_blend_flux=0.5)
        self._run_q_flux_test_arbitrary(fix_q_flux=2.1)
        self._run_q_flux_test_arbitrary(fix_q_flux=1.4, fix_blend_flux=0.25)
        self._run_q_flux_test_arbitrary(fix_q_flux=1.4, fix_blend_flux=0.)


def test_fit_fluxes():
    """
    test that when the model is updated, and fit fluxes is re-run, the fluxes
    actually change.
    """

    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1.0
    f_b = 0.5
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=False)
    #   Before update or fit_fluxes is run, chi2_per_point should be None
    assert (my_fit.chi2_per_point is None)
    my_fit.update()
    #   After update is run, chi2_per_point should have some values
    assert (len(my_fit.chi2_per_point) == 1000)
    f_s_1 = my_fit.source_flux
    chi2_1 = my_fit.chi2

    t_E_2 = pspl.parameters.t_E / (f_s + f_b)
    u_0_2 = pspl.parameters.u_0 / (f_s + f_b)
    new_model = mm.Model(
        {'t_0': pspl.parameters.t_0, 'u_0': u_0_2, 't_E': t_E_2})
    my_fit.model = new_model
    my_fit.fix_blend_flux = 0.
    my_fit.fit_fluxes()

    assert (f_s_1 != my_fit.source_flux)
    assert (chi2_1 == my_fit.chi2)

    my_fit.update()
    assert (chi2_1 != my_fit.chi2)


def create_0939_parallax_model():
    """Create Model instance with parallax"""
    model_parameters = {
        't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
        'pi_E_N': -0.248, 'pi_E_E': 0.234, 't_0_par': 2456836.2}
    coords = "17:47:12.25 -21:22:58.2"
    model_with_par = mm.Model(model_parameters, coords=coords)
    model_with_par.parallax(satellite=True, earth_orbital=True,
                            topocentric=False)
    return model_with_par


def test_satellite_and_annual_parallax_calculation():
    """
    test that data magnifications are correctly retrieved for Spitzer data.
    """
    model_with_par = create_0939_parallax_model()

    # Load Spitzer data and answers
    data_Spitzer = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True, usecols=[5])

    # Test FitData.data_magnification()
    my_fit = mm.FitData(dataset=data_Spitzer, model=model_with_par)
    ratio = my_fit.get_data_magnification() / ref_Spitzer
    almost(ratio, [1.]*len(ratio), decimal=4)


def test_get_d_A_d_params_for_point_lens_model():
    """
    Test that calculating derivatives with an ephemeris file is different from
    without an ephemeris file.
    """
    parameters = ['pi_E_N', 'pi_E_E']
    model_with_par = create_0939_parallax_model()

    data_ephm = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    fit_ephm = mm.FitData(dataset=data_ephm, model=model_with_par)
    derivs_ephm = fit_ephm.get_d_A_d_params_for_point_lens_model(parameters)

    data_no_ephm = mm.MulensData(file_name=SAMPLE_FILE_03)
    fit_no_ephm = mm.FitData(dataset=data_no_ephm, model=model_with_par)
    derivs_no_ephm = fit_no_ephm.get_d_A_d_params_for_point_lens_model(
        parameters)

    for param in parameters:
        ratio = derivs_ephm[param] / derivs_no_ephm[param]
        assert (np.abs(ratio - 1.) > 0.001).all()


def test_bad_data():
    """
    test how chi2 and chi2_per_point are affected if some datapoints are set
    to bad.

    Effectively tests
        update()
        fit_fluxes()
        get_data_magnification()
        get_model_fluxes()
        chi2
        chi2_per_point

    """

    # test that chi2 changes when using all data points vs. eliminating the
    # planet.
    (t_planet_start, t_planet_stop) = (2460982., 2460985.)
    data = mm.MulensData(file_name=SAMPLE_FILE_04_WF)
    flag_planet = (data.time > t_planet_start) & (
        data.time < t_planet_stop) | np.isnan(data.err_mag)
    data_bad = mm.MulensData(file_name=SAMPLE_FILE_04_WF, bad=flag_planet)

    (t_0, u_0, t_E) = (2460962.36458, 0.411823, 22.8092)
    point_lens_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    fit_all = mm.FitData(dataset=data, model=point_lens_model)
    fit_bad = mm.FitData(dataset=data_bad, model=point_lens_model)
    assert (fit_all.chi2 is None)
    fit_all.update()
    fit_bad.update()
    chi2_all = fit_all.chi2
    chi2_bad = fit_bad.chi2
    assert (chi2_all > chi2_bad)

    # test whether chi2_per_point is calculated for bad points.
    # not calculated --> magnification = 0, model_flux --> f_blend, dchi2=large
    # update: bad not specified --> not calculated
    # Likewise, do these tests for get_model_magnitudes
    # points:
    #   during anomaly 13055
    #   before anomaly, but excluded: 12915
    #   before anomaly, but included: 12900
    good_pt = 12900
    bad_pt = 12915
    assert (fit_bad.chi2_per_point[bad_pt] / fit_bad.chi2_per_point[good_pt] >
            100.)
    expected_mag = mm.Utils.get_mag_from_flux(fit_bad.blend_flux)
    almost(fit_bad.get_model_magnitudes()[bad_pt], expected_mag)

    # update: bad=True --> calculated
    fit_bad.update(bad=True)
    assert (fit_bad.chi2_per_point[bad_pt] / fit_bad.chi2_per_point[good_pt] <
            10.)
    almost(fit_bad.get_model_magnitudes()[bad_pt], 19.27, decimal=1)

    # update: bad=False --> not calculated
    fit_bad.update(bad=False)
    assert (fit_bad.chi2_per_point[bad_pt] / fit_bad.chi2_per_point[good_pt] >
            100.)
    almost(fit_bad.get_model_magnitudes()[bad_pt], expected_mag)

    # Test fitted fluxes are different with and without bad data points.
    assert (fit_all.source_flux > fit_bad.source_flux)


def test_scale_fluxes():
    """Specify a source_flux, blend_flux and make sure it works"""

    # Original Flux values
    f_s = 1.0
    f_b = 0.5

    # Generate fake data from a fake model
    pspl, t, A = generate_model()
    f_mod = f_s * A + f_b
    data = generate_dataset(f_mod, t)

    fit = mm.FitData(dataset=data, model=pspl)
    fit.fit_fluxes()

    num = 100
    # Test the same
    (new_flux, new_err) = fit.scale_fluxes(source_flux=f_s, blend_flux=f_b)
    almost(data.flux[num], new_flux[num])

    # Test Different
    (f_s_new, f_b_new) = (0.1, 0.)
    exp_flux = (data.flux - f_b) * f_s_new / f_s + f_b_new
    exp_err = data.err_flux * f_s_new / f_s
    (new_flux, new_err) = fit.scale_fluxes(
        source_flux=f_s_new, blend_flux=f_b_new)
    assert np.abs(data.flux[num] - new_flux[num]) > 0.5
    almost(exp_flux / new_flux, 1.)
    almost(exp_err / new_err, 1.)


class TestGetResiduals(unittest.TestCase):
    """
    test get_residuals():
    Test all keywords:
        phot_fmt: 'mag', 'flux'
        phot_fmt: 'scaled' and source_flux, blend_flux specified
        bad: True, False
    test values of residuals and errorbars
    """

    def setUp(self):
        self.model = mm.Model(
            {'t_0': 8000., 'u_0': 0.3, 't_E': 25.})
        self.generate_fake_dataset()
        self.fit = mm.FitData(model=self.model, dataset=self.dataset)
        self.fit.fit_fluxes()

    def generate_fake_dataset(self):
        """
        create a fake, perfect dataset, but with a few known outliers and
        errorbar variations.
        """
        self.dataset_properties = {
            'f_source': 10, 'f_blend': 3.5, 'errorbar': 1.}

        # Generate perfect data
        n = 3
        dt = 1.0
        times = np.arange(
            self.model.parameters.t_0 - n * self.model.parameters.t_E,
            self.model.parameters.t_0 + n * self.model.parameters.t_E,
            dt)
        flux = (self.dataset_properties['f_source'] *
                self.model.get_magnification(times) +
                self.dataset_properties['f_blend'])
        err = np.zeros(len(times)) + self.dataset_properties['errorbar']
        bad = np.zeros(len(times), dtype=bool)

        # Add outliers
        self.outliers = {'index': np.arange(0, len(times)-5, 10)+3}
        self.outliers['values'] = 10 + np.zeros(len(self.outliers['index']))
        for i in np.arange(len(self.outliers['index'])):
            if i % 5 == 0:
                self.outliers['values'][i] *= -1

            flux[self.outliers['index'][i]] += self.outliers['values'][i]
            bad[self.outliers['index'][i]] = True

        # Add errorbar variations
        self.big_errors = {'index': np.arange(0, len(times)-6, 21) + 4}
        self.big_errors['values'] = 5. + np.zeros(
            len(self.big_errors['index']))
        for i in np.arange(len(self.big_errors['index'])):
            err[self.big_errors['index'][i]] = self.big_errors['values'][i]

        assert np.sum(err) > len(err) * self.dataset_properties['errorbar']

        # Create final dataset
        self.dataset = mm.MulensData(
            [times, flux, err], phot_fmt='flux', bad=bad)

    def test_bad_keyword(self):
        """
        If bad = False, the magnification should be zero. Therefore, the flux
        calculated for the bad data points should be f_blend. If bad=True,
        the values should be the true values of the residuals.
        """
        # Bad = False
        (residuals, res_errors) = self.fit.get_residuals(
            phot_fmt='flux', bad=False)

        for index in self.outliers['index']:
            exp_residual = (self.dataset.flux[index] -
                            self.dataset_properties['f_blend'])
            almost(residuals[index], exp_residual)

        # Check errorbars
        almost(res_errors, self.dataset.err_flux)

        # Bad = True
        (residuals, res_errors) = self.fit.get_residuals(
            phot_fmt='flux', bad=True)

        for i, index in enumerate(self.outliers['index']):
            exp_residual = self.outliers['values'][i]
            almost(residuals[index], exp_residual)

        # Check errorbars
        almost(res_errors, self.dataset.err_flux)

    def test_photfmt_mag(self):
        """ check phot_fmt = 'mag' ."""
        # Bad = True
        msg = '"mag" returns residuals in the original data flux system.'
        with pytest.warns(UserWarning, match=msg):
            (residuals, res_errors) = self.fit.get_residuals(
                phot_fmt='mag', bad=True)

        # Simple sign check
        for i, index in enumerate(self.outliers['index']):
            if self.outliers['values'][i] > 0:
                assert residuals[index] < 0
            else:
                assert residuals[index] > 0

        # Value check
        for i in np.arange(len(self.dataset.time)):
            if i in self.outliers['index']:
                index = np.where(self.outliers['index'] == i)
                f_0 = self.dataset.flux[i] - self.outliers['values'][index]
                f_obs = self.dataset.flux[i]
                delta_mag = -2.5*np.log10(f_obs / f_0)
                almost(delta_mag, residuals[i])
            else:
                # Non-outliers should have zero residual
                almost(residuals[i], 0)

        # Check errorbars
        almost(res_errors, self.dataset.err_mag)

    def test_photfmt_scaled_1(self):
        """ check phot_fmt='scaled' """
        f_source_0 = 1.0
        f_blend_0 = 0.1

        # Bad = True
        (residuals, res_errors) = self.fit.get_residuals(
            phot_fmt='scaled', source_flux=f_source_0, blend_flux=f_blend_0,
            bad=True)

        model_flux = (f_source_0 *
                      self.model.get_magnification(self.dataset.time) +
                      f_blend_0)
        model_mag = mm.Utils.get_mag_from_flux(model_flux)
        for i in np.arange(len(self.dataset.time)):
            exp_flux = (f_source_0 *
                        (self.dataset.flux[i] -
                         self.dataset_properties['f_blend']) /
                        self.dataset_properties['f_source'] + f_blend_0)
            if i in self.outliers['index']:
                exp_mag = mm.Utils.get_mag_from_flux(exp_flux)
                exp_delta_mag = exp_mag - model_mag[i]
                almost(exp_delta_mag, residuals[i])
            else:
                # Non-outliers should have zero residual
                almost(residuals[i], 0)

            # Check errorbars
            exp_err_flux = (f_source_0 * self.dataset.err_flux[i] /
                            self.dataset_properties['f_source'])
            exp_err_mag = 2.5 * exp_err_flux / exp_flux / np.log(10.)
            almost(exp_err_mag, res_errors[i])
            assert self.dataset.err_mag[i] != res_errors[i]

    def test_photfmt_scaled_2(self):
        """ check phot_fmt='scaled'; true values of f_source, f_blend should
        yield errorbars identical to the true values."""
        f_source_0 = self.dataset_properties['f_source']
        f_blend_0 = self.dataset_properties['f_blend']

        # Bad = True
        (residuals, res_errors) = self.fit.get_residuals(
            phot_fmt='scaled', source_flux=f_source_0, blend_flux=f_blend_0,
            bad=True)
        almost(res_errors, self.dataset.err_mag)


class TestGradient(unittest.TestCase):
    def test_no_gradient_for_xallarap(self):
        """
        Make sure that gradient for xallarap models in not implemented.
        """
        data = mm.MulensData(file_name=SAMPLE_FILE_02)
        model = mm.Model({
            't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
            'xi_period': 100., 'xi_semimajor_axis': 0.5, 'xi_Omega_node': 90.,
            'xi_inclination': 90., 'xi_argument_of_latitude_reference': 90.})
        fit = mm.FitData(model, data)

        with self.assertRaises(NotImplementedError):
            fit.get_chi2_gradient(['t_0', 'u_0', 't_E'])


class TestFSPLGradient(unittest.TestCase):
    """ Compares various parts of the FSPL Derivative calculations to the
    results from sfit."""

    def setUp(self):
        # Read in sfit comparison file, split by dataset
        self.filenames = ['FSPL_par_Obs_1_I.pho', 'FSPL_par_Obs_2_V.pho']
        self.sfit_derivs = np.genfromtxt(
            SAMPLE_FILE_FSPL_61, dtype=None,
            names=['nob', 'k', 't', 'dAdrho', 'mag', 'db0', 'db1'])
        self._read_sfit()
        self._create_model()
        self._set_datasets()
        self._set_fits()
        self._set_indices()

    def _read_sfit(self):
        """ read in the input parameters and output matrices from sfit"""
        self.sfit_mat = FortranSFitFile51(SAMPLE_FILE_FSPL_51)
        self.sfit_mat.a[0] += 2450000.

    def _create_model(self):
        """ Initialize a model to match sfit parameters """
        parameters = ['t_0', 'u_0', 't_E', 'rho']
        self.sfit_model = mm.Model(dict(zip(parameters, self.sfit_mat.a)))
        t_star = (
            self.sfit_model.parameters.rho * self.sfit_model.parameters.t_E)
        n_t_star = 9.
        self._t_lim_1 = self.sfit_model.parameters.t_0 - n_t_star * t_star
        self._t_lim_2 = self.sfit_model.parameters.t_0 + n_t_star * t_star
        n_t_star_2 = 50.
        self._t_lim_3 = self.sfit_model.parameters.t_0 - n_t_star_2 * t_star
        self._t_lim_4 = self.sfit_model.parameters.t_0 + n_t_star_2 * t_star
        self.sfit_model.set_magnification_methods(
            [self._t_lim_3, 'finite_source_uniform_Gould94',
             self._t_lim_1, 'finite_source_LD_Yoo04',
             self._t_lim_2, 'finite_source_uniform_Gould94',
             self._t_lim_4])
        self.sfit_model.set_limb_coeff_gamma('I', self.sfit_mat.a[4])
        self.sfit_model.set_limb_coeff_gamma('V', self.sfit_mat.a[5])

    def _set_datasets(self):
        """ Read in datasets for test"""
        self.datasets = []
        for filename in self.filenames:
            bandpass = filename.split('.')[0][-1]
            mag_data = np.genfromtxt(
                join(dir_4, filename), dtype=None, encoding='utf-8',
                names=['time', 'mag', 'err'])
            (flux, err) = mm.Utils.get_flux_and_err_from_mag(
                mag_data['mag'], mag_data['err'], zeropoint=18.)
            dataset = mm.MulensData(
                [mag_data['time'], flux, err], phot_fmt='flux',
                bandpass=bandpass)
            self.datasets.append(dataset)

    def _set_fits(self):
        """ Set up fits for each individual dataset."""
        self.fits = []
        self.zs = []  # z = u / rho for each data epoch
        self.indices = []  # restrict to points affected by FS effects
        self.sfit_indices = []  # select the right parts of the sfit comparison
        # file
        for (i, dataset) in enumerate(self.datasets):
            fit = mm.FitData(
                dataset=dataset, model=self.sfit_model,
                fix_source_flux=self.sfit_mat.a[9 + i * 3],
                fix_blend_flux=self.sfit_mat.a[9 + i * 3 + 1])
            fit.fit_fluxes()
            self.fits.append(fit)

            index = ((dataset.time > self._t_lim_1) &
                     (dataset.time < self._t_lim_2))
            self.indices.append(index)

            sfit_index = np.where(self.sfit_derivs['nob'] == i + 1)
            self.sfit_indices.append(sfit_index)

            trajectory = fit.model.get_trajectory(dataset.time)
            u = np.sqrt(trajectory.x**2 + trajectory.y**2)
            z = u / self.sfit_model.parameters.rho
            self.zs.append(z)

    def _set_indices(self):
        z_break = 1.3
        zs_1_margin = 0.003
        self._indexes = []
        self._indices_not_near_1 = []
        self._indices_not_near_1_db = []
        for (zs, indices) in zip(self.zs, self.indices):
            index_large = (zs > z_break)
            index_small = (zs <= z_break)
            self._indexes.append([index_large, index_small])
            # The sfit code is not accurate near 1.0.
            near_1 = (np.abs(zs - 1.) > zs_1_margin)
            self._indices_not_near_1.append(indices & near_1)
            near_1_db = (zs < 0.88) | (zs > 1.1)
            self._indices_not_near_1_db.append(indices & near_1_db)

    def _db0_test(self, i):
        """ Test that B0prime is calculated correctly"""
        # private function check
        sfit_db0 = self.sfit_derivs[self.sfit_indices[i]]['db0']
        kwargs_ = [{'atol': 0.0005}, {'rtol': 0.01}]
        for (condition, kwargs) in zip(self._indexes[i], kwargs_):
            index_i = condition & self._indices_not_near_1_db[i]
            z = self.zs[i][index_i]
            db0 = mm.B0B1Utils().interpolate_B0prime(z)
            assert_allclose(db0, sfit_db0[index_i], **kwargs)

    def test_db0_0(self):
        """ Check that B0prime is calculated correctly for dataset 0"""
        self._db0_test(0)

    def test_db0_1(self):
        """ Check that B0prime is calculated correctly for dataset 1"""
        self._db0_test(1)

    def _db1_test(self, i):
        """ Check that B1prime is calculated correctly"""
        # private function check
        sfit_db1 = self.sfit_derivs[self.sfit_indices[i]]['db1']
        kwargs_ = [{'atol': 0.001}, {'rtol': 0.05}]
        for (condition, kwargs) in zip(self._indexes[i], kwargs_):
            index_i = condition & self._indices_not_near_1_db[i]
            z = self.zs[i][index_i]
            db1 = mm.B0B1Utils().interpolate_B1prime(z)
            assert_allclose(db1, sfit_db1[index_i], **kwargs)

    def test_db1_0(self):
        """ Check that B1prime is calculated correctly for dataset 0"""
        self._db1_test(0)

    def test_db1_1(self):
        """ Check that B1prime is calculated correctly for dataset 1"""
        self._db1_test(1)

    def _mags_test(self, i):
        """ Check that magnification is calculated correctly"""
        mags = self.fits[i].get_data_magnification()
        sfit_mags = self.sfit_derivs[self.sfit_indices[i]]['mag']
        assert_allclose(mags, sfit_mags, rtol=0.005)

    def test_mags_0(self):
        """ Check that magnification is calculated correctly for dataset 0"""
        self._mags_test(0)

    def test_mags_1(self):
        """ Check that magnification is calculated correctly for dataset 1"""
        self._mags_test(1)

    def _dA_drho_test(self, i):
        """ Check that dA / drho is calculated correctly"""
        # compare da_drho
        fs = self.fits[i].source_flux
        derivs = self.fits[i].get_d_A_d_params_for_point_lens_model(
            ['rho'])
        dA_drho = fs * derivs['rho']
        sfit_da_drho = self.sfit_derivs[self.sfit_indices[i]]['dAdrho']
        mask = self._indices_not_near_1[i]

        assert_allclose(dA_drho[mask], sfit_da_drho[mask], rtol=0.015)

    def test_dAdrho_0(self):
        """ Check that dA / drho is calculated correctly for dataset 0"""
        self._dA_drho_test(0)

    def test_dAdrho_1(self):
        """ Check that dA / drho is calculated correctly for dataset 1"""
        self._dA_drho_test(1)

    def _dA_drho_test_PLMagnification(self, i):
        """ Check that dA / drho is calculated correctly"""
        # compare da_drho
        fs = self.fits[i].source_flux
        traj = self.fits[i].get_dataset_trajectory()
        pl = mm.FiniteSourceLDYoo04Magnification(
            trajectory=traj, gamma=self.fits[i].gamma)

        sfit_da_drho = self.sfit_derivs[self.sfit_indices[i]]['dAdrho']
        mask = self._indices_not_near_1[i]

        dA_drho = fs * pl.get_d_A_d_rho()
        assert_allclose(dA_drho[mask], sfit_da_drho[mask], rtol=0.015)

        dA_drho_params = fs * pl.get_d_A_d_params(['t_0', 'rho'])['rho']
        assert_allclose(dA_drho_params[mask], sfit_da_drho[mask], rtol=0.015)

    def test_dAdrho_PLMagnification_0(self):
        """ Check that dA / drho is calculated correctly for dataset 0"""
        self._dA_drho_test_PLMagnification(0)

    def test_dAdrho_PLMagnification_1(self):
        """ Check that dA / drho is calculated correctly for dataset 1"""
        self._dA_drho_test_PLMagnification(1)

    def _set_limb_coeffs(self, model):
        for band in ['I', 'V']:
            gamma = self.sfit_model.get_limb_coeff_gamma(band)
            model.set_limb_coeff_gamma(band, gamma)

    def test_FSPL_Derivatives_tstar(self):
        """
        Make sure that FSPL Derivatives fails for models defined with tstar.
        """
        model = mm.Model(
            {'t_0': self.sfit_model.parameters.t_0,
             'u_0': self.sfit_model.parameters.u_0,
             't_E': self.sfit_model.parameters.t_E,
             't_star': self.sfit_model.parameters.t_star})
        self._set_limb_coeffs(model)
        fit = mm.FitData(model=model, dataset=self.datasets[0])

        with self.assertRaises(AttributeError):
            fit.get_d_A_d_rho()

    def test_check_FSPLDerivs_errors_1(self):
        parameters = ['t_0', 'u_0', 't_E', 'rho']
        model = mm.Model(dict(zip(parameters, self.sfit_mat.a)))
        self._set_limb_coeffs(model)

        t_star = model.parameters.rho * model.parameters.t_E
        n_t_star = 9.
        t_lim_1 = model.parameters.t_0 - n_t_star * t_star
        t_lim_2 = model.parameters.t_0 + n_t_star * t_star
        model.set_magnification_methods([t_lim_1, 'finite_source_LD_WittMao94', t_lim_2])
        fit = mm.FitData(model=model, dataset=self.datasets[0])
        with self.assertRaises(NotImplementedError):
            fit.get_d_A_d_rho()

    def test_magnification_methods_parameters(self):
        parameters = ['t_0', 'u_0', 't_E', 'rho']
        model = mm.Model(dict(zip(parameters, self.sfit_mat.a)))
        t_star = model.parameters.rho * model.parameters.t_E
        n_t_star = 9.
        t_lim_1 = model.parameters.t_0 - n_t_star * t_star
        t_lim_2 = model.parameters.t_0 + n_t_star * t_star
        model.set_magnification_methods(
            [t_lim_1, 'finite_source_uniform_Gould94', t_lim_2])
        with self.assertRaises(KeyError):
            model.set_magnification_methods_parameters(
                {'vbbl': {'accuracy': 0.005}})

        model.set_magnification_methods_parameters(
            {'finite_source_uniform_Gould94': {'accuracy': 0.005}})
        with self.assertRaises(ValueError):
            model.get_magnification(np.arange(t_lim_1, t_lim_2, 0.1))


class TestFSPLGradient2(TestFSPLGradient):

    def setUp(self):
        fort_61 = join(dir_4, 'test_2', 'fort.61')
        fort_62 = join(dir_4, 'test_2', 'fort.62')
        fort_51 = join(dir_4, 'test_2', 'fort.51')
        self.filenames = [join('test_2', 'FSPL_Obs_1_I.pho'),
                          join('test_2', 'FSPL_Obs_2_V.pho')]

        self.sfit_derivs = np.genfromtxt(
            fort_61, dtype=None, encoding='utf-8',
            names=['nob', 'k', 't', 'dAdrho', 'mag', 'db0', 'db1'])
        self.sfit_partials = np.genfromtxt(
            fort_62, dtype=None, encoding='utf-8',
            names=['nob', 'k', 't', 'dfdt0', 'dfdu0', 'dfdtE', 'dfdrho',
                   'dAdu', 'df', 'res', 'sig2'])
        self.sfit_mat = FortranSFitFile51(fort_51)
        self.sfit_mat.a[0] += 2450000.

        self._create_model()
        self._set_datasets()
        self._set_fits()
        self._set_indices()

    def test_dA_dparams(self):
        params = ['t_0', 'u_0', 't_E', 'rho']
        for i, fit in enumerate(self.fits):
            dA_dparam = fit.get_d_A_d_params_for_point_lens_model(params)
            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                df_dparam = fit.source_flux * dA_dparam[param]
                sfit_df_dparam = self.sfit_partials[
                    self.sfit_indices[i]]['dfd{0}'.format(short_param)]
                mask = self._indices_not_near_1[i]
                assert_allclose(
                    df_dparam[mask], sfit_df_dparam[mask], rtol=0.015)

    def test_chi2_gradient(self):
        params = ['t_0', 'u_0', 't_E', 'rho']
        for i, fit in enumerate(self.fits):
            gradient = fit.get_chi2_gradient(params)
            sfit_res = self.sfit_partials[self.sfit_indices[i]]['res']

            # Test residuals from model & Errors
            res, err = fit.get_residuals(phot_fmt='flux')
            index_peak = (self.zs[i] < 15.)
            index_wings = (self.zs[i] > 15.)
            assert_allclose(res[index_peak], sfit_res[index_peak], rtol=0.01)
            assert_allclose(res[index_wings], sfit_res[index_wings], atol=0.01)
            sig2 = self.sfit_partials[self.sfit_indices[i]]['sig2']
            assert_allclose(err**2, sig2, rtol=0.01)

            # Test gradient
            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                partials = self.sfit_partials[
                    self.sfit_indices[i]]['dfd{0}'.format(short_param)]
                sfit_gradient = np.sum(-2. * sfit_res * partials / sig2)
                assert_allclose(gradient[j], sfit_gradient, rtol=0.01)

    def test_dAdu_FSPLError(self):
        for i, fit in enumerate(self.fits):
            with self.assertRaises(NotImplementedError):
                fit.get_d_A_d_u_for_PSPL_model()
            with self.assertRaises(NotImplementedError):
                fit.get_d_A_d_u_for_FSPL_model()
            with self.assertRaises(NotImplementedError):
                fit.get_d_A_d_u_for_point_lens_model()

    def test_dAdu_PSPL(self):
        params = ['t_0', 'u_0', 't_E']
        sfit_PSPL_model = mm.Model(dict(zip(params, self.sfit_mat.a)))
        for (i, dataset) in enumerate(self.datasets):
            fit = mm.FitData(
                dataset=dataset, model=sfit_PSPL_model,
                fix_source_flux=self.sfit_mat.a[9 + i * 3],
                fix_blend_flux=self.sfit_mat.a[9 + i * 3 + 1])
            fit.fit_fluxes()

            pl = mm.PointSourcePointLensMagnification(
                trajectory=fit.get_dataset_trajectory())
            dAdu = pl.get_d_A_d_u()
            assert_allclose(dAdu, self.sfit_partials[self.sfit_indices[i]]['dAdu'], rtol=0.005)


def test_FSPLDerivs_get_satellite_coords():
    """Test that satellite_skycoord propagate correctly through the code."""
    # Inputs
    times = [2456445.0, 2457328.0]
    mags = [16., 15.]
    errs = [0.01, 0.02]

    # Expected results
    ra_1 = 15 * (8 + 26 / 60. + 37.19 / 3600.)
    dec_1 = 18 + 30 / 60. + 37.4 / 3600.
    ra_2 = 15 * (17 + 40 / 60. + 4.98 / 3600.)
    dec_2 = -23 - 26 / 60. - 38.2 / 3600.

    dataset = mm.MulensData(
        [times, mags, errs], phot_fmt='mag',
        ephemerides_file=SAMPLE_FILE_03_EPH)
    model = mm.Model({'t_0': 2457000., 'u_0': 0.01, 't_E': 100., 'rho': 0.02})

    # Partial test with _direct method
    model.default_magnification_method = 'finite_source_uniform_Gould94_direct'
    fit = mm.FitData(dataset=dataset, model=model)
    fit._set_data_magnification_curves()
    mag_curve = fit._data_magnification_curve
    mag_curve._set_point_lens_magnification_objects()
    assert len(mag_curve._magnification_objects) == 1

    # Continue original test without the _direct method
    model.default_magnification_method = 'finite_source_uniform_Gould94'
    fit = mm.FitData(dataset=dataset, model=model)
    fit._set_data_magnification_curves()
    mag_curve = fit._data_magnification_curve
    mag_curve._set_point_lens_magnification_objects()
    assert len(mag_curve._magnification_objects) == 1

    for derivs_obj in mag_curve._magnification_objects.values():
        if times[0] in derivs_obj.trajectory.times:
            result_1 = derivs_obj.trajectory.satellite_skycoord[0]
            np.testing.assert_almost_equal(result_1.ra.value, ra_1, decimal=3)
            np.testing.assert_almost_equal(result_1.dec.value, dec_1, decimal=3)

        if times[-1] in derivs_obj.trajectory.times:
            result_2 = derivs_obj.trajectory.satellite_skycoord[-1]
            np.testing.assert_almost_equal(result_2.ra.value, ra_2, decimal=3)
            np.testing.assert_almost_equal(result_2.dec.value, dec_2, decimal=3)


def test_get_trajectory_1L2S_satellite_parallax():
    """test parallax calculation with Spitzer data"""
    model_parameters = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
                        'pi_E_N': -0.248, 'pi_E_E': 0.234,
                        't_0_par': 2456836.2}
    coords = "17:47:12.25 -21:22:58.2"

    model_with_par = mm.Model(model_parameters, coords=coords)
    model_with_par.parallax(satellite=True, earth_orbital=True,
                            topocentric=False)

    data_Spitzer = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)

    fit = mm.FitData(model=model_with_par, dataset=data_Spitzer)

    ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True)

    trajectory = fit.get_dataset_trajectory()

    ratio_x = trajectory.x / ref_Spitzer[6]
    ratio_y = trajectory.y / ref_Spitzer[7]
    np.testing.assert_almost_equal(ratio_x, 1., decimal=2)
    np.testing.assert_almost_equal(ratio_y, 1., decimal=3)


def test_bad_data_w_ephemerides():
    """
    Test that satellite_skycoords is correctly masked when there are bad data.
    """
    model_with_par = create_0939_parallax_model()

    # Load Spitzer data and answers
    data_Spitzer = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    bad = np.zeros(len(data_Spitzer.time), dtype=bool)
    bad[-10:] = True
    data_Spitzer.bad = bad

    my_fit = mm.FitData(dataset=data_Spitzer, model=model_with_par)
    magnification = my_fit.get_data_magnification(bad=False)
    np.testing.assert_almost_equal(magnification[-10:], 0.)


class Test2L2S(unittest.TestCase):

    def setUp(self):
        self.model_type = '2L2S'
        self._setup()

    def _setup(self):
        filename = 'OB151489_simulated_data_{0}.dat'.format(self.model_type)
        file_path = join(dir_2, filename)
        self.header_info = self.parse_header(file_path)

        self.data = mm.MulensData(
            file_name=file_path, phot_fmt=self.header_info['phot_fmt'],)

        self.model = mm.Model(self.header_info['params'])
        for param in self.header_info['params']:
            if 'rho' in param:
                num = param.split('_')[-1]
                self.model.set_magnification_methods(self.header_info['mag_methods'], int(num))

    def parse_header(self, file_path):
        header_info = {}
        params = {}
        param_heads = ['t_0', 'u_0', 'rho']
        with open(file_path) as file_:
            for line in file_.readlines():
                if '# Source ' in line:
                    n = line[9:10]
                    source_params = ast.literal_eval(line[11:].strip())
                    for key, value in source_params.items():
                        if key in param_heads:
                            params['{0}_{1}'.format(key, n)] = value
                        elif key not in params.keys():
                            params[key] = value
                elif '# source:' in line:
                    fluxes = line[9:].split()
                    header_info['source fluxes'] = [float(flux) for flux in fluxes]
                elif '# blend:' in line:
                    header_info['blend flux'] = float(line[8:])
                elif '# mag_methods:' in line:
                    header_info['mag_methods'] = ast.literal_eval(line[15:].strip())
                elif 'Date,' in line:
                    header_info['phot_fmt'] = line.split(',')[1].strip().lower()

        header_info['params'] = params
        return header_info

    def _test_fitted_fluxes(self, fit):
        fit.fit_fluxes()
        np.testing.assert_almost_equal(
            fit.blend_flux, self.header_info['blend flux'], decimal=3)
        for i in range(self.model.n_sources):
            np.testing.assert_almost_equal(
                fit.source_fluxes[i], self.header_info['source fluxes'][i],
                decimal=3)

    def test_multiple_source(self):
        """
        Test that we can fit fluxes for N sources.
        """
        self._test_fitted_fluxes(mm.FitData(dataset=self.data, model=self.model))

    def test_fix_source_flux(self):
        # Each source separately
        for i, flux in enumerate(self.header_info['source fluxes']):
            fix_source_flux = [False for i in range(len(self.header_info['source fluxes']))]
            fix_source_flux[i] = flux
            self._test_fitted_fluxes(mm.FitData(dataset=self.data, model=self.model, fix_source_flux=fix_source_flux))

        # All sources
        self._test_fitted_fluxes(
            mm.FitData(dataset=self.data, model=self.model, fix_source_flux=self.header_info['source fluxes']))

    def test_fix_source_flux_arbitrary(self):
        fix_source_flux = [0.04, 0.01, 0.005]
        fit = mm.FitData(
            dataset=self.data, model=self.model,
            fix_source_flux=fix_source_flux[0:self.model.n_sources])

        fit.fit_fluxes()
        for i in range(self.model.n_sources):
            np.testing.assert_almost_equal(
                fit.source_fluxes[i], fix_source_flux[i])

    def test_fix_source_flux_bad(self):
        fix_source_flux = [0.04, 0.01, 0.005, 16., 27.]
        with self.assertRaises(ValueError):
            mm.FitData(dataset=self.data, model=self.model, fix_source_flux=fix_source_flux)

    def test_fix_source_flux_ratio(self):
        source_flux_ratio = self.header_info['source fluxes'][1] / self.header_info['source fluxes'][0]
        fit = mm.FitData(dataset=self.data, model=self.model, fix_source_flux_ratio=source_flux_ratio)
        self._test_fitted_fluxes(fit)
        np.testing.assert_almost_equal(fit.source_flux_ratio, [source_flux_ratio])

    def test_fix_source_flux_ratio_bad(self):
        fix_source_flux_ratio = [0.04, 0.01, 0.005, 16., 27.]
        with self.assertRaises(ValueError):
            mm.FitData(dataset=self.data, model=self.model, fix_source_flux=fix_source_flux_ratio)

    def test_fix_source_flux_ratio_arbitrary(self):
        source_flux_ratios = [0.01, 0.005]
        fit = mm.FitData(
            dataset=self.data, model=self.model, fix_source_flux_ratio=source_flux_ratios[0:self.model.n_sources-1])
        fit.fit_fluxes()
        for i in range(1, self.model.n_sources):
            flux_ratio = fit.source_fluxes[i] / fit.source_fluxes[0]
            np.testing.assert_almost_equal(
                flux_ratio, source_flux_ratios[i - 1], decimal=3)
            np.testing.assert_almost_equal(fit.source_flux_ratio[i-1], flux_ratio)

    def test_fix_blend_flux(self):
        # Might as well...
        blend_flux = self.header_info['blend flux']
        self._test_fitted_fluxes(mm.FitData(dataset=self.data, model=self.model, fix_blend_flux=blend_flux))


class Test1L3S(Test2L2S):

    def setUp(self):
        self.model_type = '1L3S'
        self._setup()

    def test_fix_source_flux_ratio(self):
        source_flux_ratios = [
            self.header_info['source fluxes'][i] / self.header_info['source fluxes'][0]
            for i in range(1, self.model.n_sources)]

        for i, flux in enumerate(source_flux_ratios):
            fix_source_flux_ratio = [False for i in range(len(source_flux_ratios))]
            fix_source_flux_ratio[i] = flux
            fit = mm.FitData(dataset=self.data, model=self.model, fix_source_flux_ratio=fix_source_flux_ratio)
            self._test_fitted_fluxes(fit)
            np.testing.assert_almost_equal(fit.source_flux_ratio[i], flux)

        fit = mm.FitData(dataset=self.data, model=self.model, fix_source_flux_ratio=source_flux_ratios)
        self._test_fitted_fluxes(fit)
        np.testing.assert_almost_equal(fit.source_flux_ratio, source_flux_ratios)

# Tests to add:
#
# test get_chi2_gradient(), chi2_gradient:
#   Effectively covered by unit tests in event.py
#
# properties:
#   chi2, chi2_per_point, source_flux, source_fluxes, blend_flux, q_flux,
#   dataset, model
