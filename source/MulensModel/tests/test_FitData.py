import numpy as np
from numpy.testing import assert_almost_equal as almost
import unittest
import os.path

import MulensModel as mm

dir_1 = os.path.join(mm.DATA_PATH, 'photometry_files', 'OB140939')
dir_2 = os.path.join(mm.DATA_PATH, 'unit_test_files')
dir_3 = os.path.join(mm.DATA_PATH, 'ephemeris_files')

SAMPLE_FILE_02 = os.path.join(dir_1, 'ob140939_OGLE.dat')  # HJD'
SAMPLE_FILE_02_REF = os.path.join(dir_2, 'ob140939_OGLE_ref_v1.dat')  # HJD'
SAMPLE_FILE_03 = os.path.join(dir_1, 'ob140939_Spitzer.dat')  # HJD'
SAMPLE_FILE_03_EPH = os.path.join(dir_3, 'Spitzer_ephemeris_01.dat')  # UTC
SAMPLE_FILE_03_REF = os.path.join(dir_2, 'ob140939_Spitzer_ref_v1.dat')  # HJD'
SAMPLE_FILE_04_WF = os.path.join(mm.DATA_PATH, 'WFIRST_1827.dat')

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


def execute_test_blend_fixed(f_b):
    # test for when source flux is to be determined, but blend flux is a
    # fixed value

    pspl, t, A = generate_model()

    # secret source flux
    f_s = 1.0
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=f_b,
        fix_source_flux=False)
    my_fit.fit_fluxes()

    almost(my_fit.source_flux, f_s)


class BinarySourceTest():

    def __init__(self):
        self.f_s_1 = 1
        self.f_s_2 = 1.2
        self.f_b = 0.5

        self.setup_dataset()

    def setup_dataset(self):
        self.model, self.t, self.A_1, self.A_2 = generate_binary_model()
        f_mod = self.f_s_1 * self.A_1 + self.f_s_2 * self.A_2 + self.f_b
        self.dataset = generate_dataset(f_mod, self.t)

    def run_test(
            self, fix_blend_flux=False, fix_source_flux=False,
            fix_q_flux=False):

        self.my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_blend_flux=fix_blend_flux,
            fix_source_flux=fix_source_flux, fix_source_flux_ratio=fix_q_flux)
        self.my_fit.fit_fluxes()

        almost(self.my_fit.blend_flux, self.f_b)
        almost(self.my_fit.source_fluxes[0], self.f_s_1)
        almost(self.my_fit.source_fluxes[1], self.f_s_2)

        # Test get_model_fluxes() for 2 sources
        peak_index = 500
        mod_fluxes = self.my_fit.get_model_fluxes()
        almost(mod_fluxes[peak_index], self.dataset.flux[peak_index])


def execute_test_binary_source(q_flux=False):
    # test for when blend flux and source flux are to be determined for binary
    # sources with q-flux

    test = BinarySourceTest()
    if q_flux:
        fix_q_flux = test.f_s_2 / test.f_s_1
    else:
        fix_q_flux = False

    test.run_test(fix_q_flux=fix_q_flux)


# *** Actual tests below ***
def test_default():
    """
    test for when blend flux and source flux are to be determined
    """
    pspl, t, A = generate_model()

    # secrets
    f_s = 1.0
    f_b = 0.5
    # generate f_mod
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=False)
    my_fit.fit_fluxes()

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_flux, f_s)

    # Test get_model_fluxes() for 1 source
    peak_index = 500
    mod_fluxes = my_fit.get_model_fluxes()
    almost(mod_fluxes[peak_index], my_dataset.flux[peak_index])


def test_blend_zero():
    """
    test for when source flux is to be determined, but blend flux is zero
    """
    execute_test_blend_fixed(f_b=0.)


def test_blend_fixed():
    """
    test for when source flux is to be determined, but blend flux is
    zero
    """
    execute_test_blend_fixed(f_b=0.5)
    execute_test_blend_fixed(f_b=-0.5)


def test_source_fixed():
    """
    test for when blend flux is to be determined, but source flux is a fixed
    value
    """

    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1.0
    f_b = 0.5
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=f_s)
    my_fit.fit_fluxes()

    almost(my_fit.blend_flux, f_b)


def test_both_fixed():
    """
    test for when both fluxes are fixed --> evaluate chi2 but not fluxes.
    """
    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1.0
    f_b = 0.5
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=f_b,
        fix_source_flux=f_s)
    my_fit.update()

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_flux, f_s)


def test_binary_source():
    """Test a binary source model with all free parameters."""
    execute_test_binary_source(q_flux=False)


def test_binary_source_fixed():
    """
    Test the three cases for fixing each of the three components of a binary
    source model
    """
    test = BinarySourceTest()
    test.run_test(fix_source_flux=[False, False])
    test.run_test(fix_source_flux=[1.0, False])
    test.run_test(fix_source_flux=[False, 1.2])
    test.run_test(fix_source_flux=[1.0, 1.2])
    test.run_test(fix_blend_flux=0.5)
    test.run_test(fix_source_flux=[1.0, 1.2], fix_blend_flux=0.5)
    test.run_test(fix_source_flux=[1.0, False], fix_blend_flux=0.5)


class TestFitData(unittest.TestCase):
    def test_init_1(self):
        with self.assertRaises(ValueError):
            test = BinarySourceTest()
            test.run_test(fix_source_flux=1.0)


def test_binary_qflux():
    """
    test for when blend flux and source flux are to be determined for binary
    sources with q-flux
    """

    execute_test_binary_source(q_flux=True)


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
    assert(my_fit.chi2_per_point is None)
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

    assert(f_s_1 != my_fit.source_flux)
    assert(chi2_1 == my_fit.chi2)

    my_fit.update()
    assert(chi2_1 != my_fit.chi2)


def test_chi2_per_point():
    """Test that the chi2 shape is correct for multiple sources, i.e. = number
    of epochs, rather than epochs * sources."""
    test_object = BinarySourceTest()
    my_fit = mm.FitData(model=test_object.model, dataset=test_object.dataset)
    my_fit.update()

    assert(my_fit.chi2_per_point.shape == (test_object.dataset.n_epochs,))


def test_satellite_and_annual_parallax_calculation():
    """
    test that data magnifications are correctly retrieved for Spitzer data.
    """

    # Create Model
    model_parameters = {
        't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
        'pi_E_N': -0.248, 'pi_E_E': 0.234, 't_0_par': 2456836.2}
    coords = "17:47:12.25 -21:22:58.2"
    model_with_par = mm.Model(model_parameters, coords=coords)
    model_with_par.parallax(satellite=True, earth_orbital=True,
                            topocentric=False)

    # Load Spitzer data and answers
    data_Spitzer = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True, usecols=[5])

    # Test FitData.data_magnification()
    my_fit = mm.FitData(dataset=data_Spitzer, model=model_with_par)
    ratio = my_fit.get_data_magnification() / ref_Spitzer
    np.testing.assert_almost_equal(ratio, [1.]*len(ratio), decimal=4)


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
    assert(fit_all.chi2 is None)
    fit_all.update()
    fit_bad.update()
    chi2_all = fit_all.chi2
    chi2_bad = fit_bad.chi2
    assert(chi2_all > chi2_bad)

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

class FortranSFitFile(object):
    """
    Class to parse the comparison file generated by the Fortran version of sfit.
    """

    def __init__(self, filename):
        input_file = open(filename, 'r')
        attr = None
        for line in input_file.readlines():
            str_vec = line.strip().split()
            if str_vec[0] == '#':
                attr = str_vec[-1]
            else:
                if len(str_vec) == 1:
                    value = float(str_vec[0])
                else:
                    value = np.array([float(item) for item in str_vec])

                self.__setattr__(attr, value)

        input_file.close()


class TestFSPLGradient(unittest.TestCase):

    def setUp(self):
        # Read in sfit comparison file, split by dataset
        fspl_dir = 'fspl_derivs'
        self.sfit_derivs = np.genfromtxt(
            os.path.join(dir_2, fspl_dir, 'fort.61'), dtype=None,
            names=['nob', 'k', 't', 'dAdrho', 'mag', 'db0', 'db1'])

        # Create the model
        self.sfit_mat = FortranSFitFile(os.path.join(dir_2, fspl_dir, 'fort.51'))
        self.sfit_model = mm.Model({
            't_0': self.sfit_mat.a[0] + 2450000., 'u_0': self.sfit_mat.a[1],
            't_E': self.sfit_mat.a[2], 'rho': self.sfit_mat.a[3]})
        self.sfit_model.parameters.t_0_par = self.sfit_mat.a[0] + 2450000.
        t_star = self.sfit_model.parameters.rho * self.sfit_model.parameters.t_E
        n_t_star = 9.
        self.sfit_model.set_magnification_methods([
            self.sfit_model.parameters.t_0 - n_t_star * t_star,
            'finite_source_LD_Yoo04_direct',
            self.sfit_model.parameters.t_0 + n_t_star * t_star])
        self.sfit_model.set_limb_coeff_gamma('I', 0.44)
        self.sfit_model.set_limb_coeff_gamma('V', 0.72)

        self.filenames = ['FSPL_par_Obs_1_I.pho', 'FSPL_par_Obs_2_V.pho']
        self.datasets = []
        for filename in self.filenames:
            bandpass = filename.split('.')[0][-1]
            dataset = mm.MulensData(
                file_name=os.path.join(dir_2, fspl_dir, filename),
                phot_fmt='mag', bandpass=bandpass)
            self.datasets.append(dataset)

        self.fits = []
        self.zs = []
        self.indices = []
        self.sfit_indices = []
        for i, dataset in enumerate(self.datasets):
            fit = mm.FitData(
                dataset=dataset, model=self.sfit_model,
                fix_source_flux=self.sfit_mat.a[9 + i * 3],
                fix_blend_flux=self.sfit_mat.a[9 + i * 3 + 1])
            fit.fit_fluxes()
            self.fits.append(fit)

            index = ((dataset.time >
                      self.sfit_model.parameters.t_0 - n_t_star * t_star) &
                     (dataset.time <
                      self.sfit_model.parameters.t_0 + n_t_star * t_star))
            self.indices.append(index)

            sfit_index = np.where(self.sfit_derivs['nob'] == i + 1)
            self.sfit_indices.append(sfit_index)

            traj = fit.model.get_trajectory(dataset.time[index])
            z = (np.sqrt(traj.x ** 2 + traj.y ** 2) /
                 self.sfit_model.parameters.rho)
            self.zs = z

    def _db0_test(self, i):
        sfit_db0 = self.sfit_derivs[self.sfit_indices[i]]['db0']
        db0 = self.fits[i]._get_B0_prime(self.zs[i])
        np.testing.assert_allclose(
            db0, sfit_db0[self.indices[i]], atol=0.001)

    def test_db0_0(self):
        self._db0_test(0)

    def test_db0_1(self):
        self._db0_test(1)

    def _db1_test(self, i):
        sfit_db1 = self.sfit_derivs[self.sfit_indices[i]]['db1']
        db1 = self.fits[i]._get_B1_prime(self.zs[i])
        np.testing.assert_allclose(
            db1, sfit_db1[self.indices[i]], atol=0.001)

    def test_db1_0(self):
        self._db1_test(0)

    def test_db1_1(self):
        self._db1_test(1)

    def _mags_test(self, i):
        mags = self.fits[i].get_data_magnification()
        sfit_mags = self.sfit_derivs[self.sfit_indices[i]]['mag']
        print('shapes: z, mm mags, sfit mags')
        print(self.zs[i].shape, mags.shape, sfit_mags.shape)
        test_arr = np.vstack((
            self.zs[i], mags[self.indices[i]], sfit_mags[self.sfit_indices[i]],
            mags[self.indices[i]] / sfit_mags[self.sfit_indices[i]],
        ))
        print('z, mm mags, sfit mags, ratio')
        print(test_arr.transpose())
        np.testing.assert_allclose(
            mags, sfit_mags[self.sfit_indices[i]]['mag'], rtol=0.005)

    def test_mags_0(self):
        self._mags_test(0)

    def test_mags_1(self):
        self._mags_test(1)

    def _dA_drho_test(self, i):
        fs = self.fits[i].source_flux
        derivs = fs * self.fits[i].get_d_A_d_rho()
        # compare da_drho
        sfit_da_drho = self.sfit_derivs[self.sfit_indices[i]]['dAdrho']
        test_arr = np.vstack((
            self.zs[i], derivs[self.indices[i]], sfit_da_drho[self.indices[i]],
            derivs[self.indices[i]] / sfit_da_drho[self.indices[i]],
        ))
        print('z, mm deriv, sfit deriv, ratio')
        print(test_arr.transpose())
        np.testing.assert_allclose(
            derivs[self.indices[i]], sfit_da_drho[self.indices[i]], rtol=0.01)

    def test_dAdrho_0(self):
        self._dA_drho_test(0)

    def test_dAdrho_1(self):
        self._dA_drho_test(1)


# Tests to add:
#
# test get_chi2_gradient(), chi2_gradient:
#   Effectively covered by unit tests in event.py
#
# properties:
#   chi2, chi2_per_point, source_flux, source_fluxes, blend_flux, q_flux,
#   dataset, model
