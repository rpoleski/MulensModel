import unittest

import numpy as np
import os

import MulensModel as mm
import fortran_files
from test_FitData import create_0939_parallax_model, SAMPLE_FILE_03, SAMPLE_FILE_03_EPH

SAMPLE_FILE = os.path.join(mm.DATA_PATH, 'unit_test_files', 'FSPL_test_1.dat')
PSPL_SAMPLE_DIR = os.path.join(mm.DATA_PATH, 'unit_test_files', 'fspl_derivs', 'test_PointLensClasses')


def get_file_params(filename):
    """Read in the model parameters used to create the file"""
    with open(filename) as data_file:
        lines = data_file.readlines()
        ulens_params = lines[2].split()
    model = mm.ModelParameters({'t_0': float(ulens_params[1]), 'u_0': float(ulens_params[2]),
                                't_E': float(ulens_params[3]), 'rho': float(ulens_params[4])})
    return (model, float(ulens_params[5]))


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
    test_b_1 = test_FSPL_LD._B_1_function()
    np.testing.assert_almost_equal(test_b_1, data['b_1'], decimal=4)


def test_B_1_function_direct_and_interpolated():
    """
    FiniteSourceLDYoo04Magnification was failing at some point with normalized trajectory positions
    sqrt(x**2+y**2)/rho both smaller and larger than 10^4 (limit of interpolation).
    Here we compare magnification to PSPL with small correction
    """
    model_1 = mm.ModelParameters({'t_0': 0, 'u_0': 0.9, 't_E': 1., 'rho': 0.0001})
    model_2 = mm.ModelParameters({'t_0': 0, 'u_0': 0.9, 't_E': 1.})
    trajectory_1 = mm.Trajectory([0., 1., 1.8], model_1)
    trajectory_2 = mm.Trajectory([0., 1., 1.8], model_2)
    test_FSPL_LD = mm.FiniteSourceLDYoo04Magnification(trajectory=trajectory_1, gamma=0.44)
    test_PSPL = mm.PointSourcePointLensMagnification(trajectory=trajectory_2)
    mag_1 = test_FSPL_LD.get_magnification()[1:]
    mag_2 = test_PSPL.get_magnification()[1:]
    np.testing.assert_almost_equal(mag_1/mag_2, 1., decimal=6)


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

    mag_curve = mm.FiniteSourceUniformGould94Magnification(
        trajectory=trajectory)
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


# Make sure every element of the PointLensMagnification classes are tested.
class TestPointSourcePointLensMagnification(unittest.TestCase):

    def setUp(self):
        self.sfit_files = fortran_files.read_sfit_files(PSPL_SAMPLE_DIR)

        parameters = ['t_0', 'u_0', 't_E', 'rho']
        self.parameters = mm.ModelParameters(
            dict(zip(parameters, self.sfit_files['51'].a)))

        self.gammas = self.sfit_files['51'].a[4:5]  # Cludgy and inflexible.
        self.trajectories = []
        self._set_trajectories()

        self.mag_objs = []
        self._set_mag_objs()

    def _set_trajectories(self):
        for nob_indices in self.sfit_files['63'].sfit_nob_indices:
            trajectory = mm.Trajectory(
                self.sfit_files['63'].t[nob_indices], self.parameters)
            self.trajectories.append(trajectory)

    def _set_mag_objs(self):
        for trajectory in self.trajectories:
            mag_obj = mm.PointSourcePointLensMagnification(trajectory)
            self.mag_objs.append(mag_obj)

    def test_get_pspl_magnification(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            pspl_mag = mag_obj.get_pspl_magnification()
            np.testing.assert_allclose(
                pspl_mag, self.sfit_files['63'].amp[nob_indices], rtol=0.0001)

    def test_get_magnification(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            mag = mag_obj.get_magnification()
            np.testing.assert_allclose(
                mag, self.sfit_files['63'].amp[nob_indices], rtol=0.0001)

    def _get_factor_b1(self, nob_indices, gamma):
        fspl_factor = self.sfit_files['63'].amp[nob_indices] * (
                self.sfit_files['61'].db0[nob_indices] -
                gamma * self.sfit_files['61'].db1[nob_indices])
        fspl_factor /= self.sfit_files['51'].a[3]  # rho
        fspl_factor += self.sfit_files['62'].dAdu[nob_indices] * (
                self.sfit_files['63'].b0[nob_indices] -
                gamma * self.sfit_files['63'].b1[nob_indices])

        return fspl_factor

    def test_get_d_A_d_params(self):
        """
        df/dparams = fs * dA/dparams (FSPL)

        dA/dparams (PSPL) = d_A_d_u * d_u_d_params[key]
        dA/dparams (FSPL) = d_u_d_params[key] * factor =
            factor * dA/dparams(PSPL) / dA_du

        dA_dparams(PSPL) = df/dparams * dA_du / fs / factor
        """
        params = ['t_0', 'u_0', 't_E']

        for (nob_indices, source_flux, gamma, mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.gammas, self.mag_objs):

            dA_dparam = mag_obj.get_d_A_d_params(params)

            fspl_factor = self._get_factor_b1(nob_indices, gamma)

            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                sfit_df_dparam = self.sfit_files['62'].data[
                    'dfd{0}'.format(short_param)][nob_indices]
                sfit_dA_dparam = (sfit_df_dparam *
                                  self.sfit_files['62'].dAdu[nob_indices] /
                                  source_flux / fspl_factor)
                np.testing.assert_allclose(
                    dA_dparam[param], sfit_dA_dparam, rtol=0.015)

    def test_get_d_u_d_params(self):
        """
        PSPL
        d_A_d_params[key] = d_A_d_u * d_u_d_params[key]

        FSPL
        if key == 'rho':
            d_A_d_params[key] = self.get_d_A_d_rho()
        else:
            d_A_d_params[key] = d_u_d_params[key] * factor

        sfit returns: FSPL:
            61 dA/drho
            62 df/dparams, dAdu

        df/dparams = fs * dA/dparams (FSPL)

        So, df/dparams / fs / factor = d_u_d_params
        """

        params = ['t_0', 'u_0', 't_E']

        for (nob_indices, source_flux, gamma, mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.gammas, self.mag_objs):

            du_dparam = mag_obj.get_d_u_d_params(params)

            fspl_factor = self._get_factor_b1(nob_indices, gamma)

            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                sfit_df_dparam = self.sfit_files['62'].data[
                    'dfd{0}'.format(short_param)][nob_indices]
                sfit_dA_dparam = sfit_df_dparam / source_flux / fspl_factor
                np.testing.assert_allclose(
                    du_dparam[param], sfit_dA_dparam, rtol=0.015)

    def test_get_d_A_d_u(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices, self.mag_objs):

            dA_du = mag_obj.get_d_A_d_u()
            np.testing.assert_allclose(
                dA_du, self.sfit_files['62'].dAdu[nob_indices], rtol=0.015)

    def test_pspl_magnification(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            np.testing.assert_allclose(
                mag_obj.pspl_magnification,
                self.sfit_files['63'].amp[nob_indices],
                rtol=0.0001)

    def test_magnification(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            with self.assertRaises(AttributeError):
                _ = mag_obj.magnification

            mag_obj.get_magnification()
            np.testing.assert_allclose(
                mag_obj.magnification,
                self.sfit_files['63'].amp[nob_indices],
                rtol=0.0001)

    def test_u_(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            np.testing.assert_allclose(
                mag_obj.u_, self.sfit_files['63'].x[nob_indices],
                rtol=0.0001)

    def test_u_2(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            np.testing.assert_allclose(
                mag_obj.u_2, self.sfit_files['63'].x2[nob_indices],
                rtol=0.0001)


class TestFiniteSourceUniformGould94Magnification(TestPointSourcePointLensMagnification):

    def setUp(self):
        TestPointSourcePointLensMagnification.setUp(self)

        self.zs = []
        for mag_obj in self.mag_objs:
            z = mag_obj.u_ / mag_obj.trajectory.parameters.rho
            self.zs.append(z)

        self._indexes = []
        self._indices_not_near_1 = []
        self._indices_not_near_1_db = []
        self.indices_mag_test = []

        self._set_indices()

    def _set_mag_objs(self):
        for trajectory in self.trajectories:
            mag_obj = mm.FiniteSourceUniformGould94Magnification(
                trajectory=trajectory)
            self.mag_objs.append(mag_obj)

    def _set_indices(self):
        z_break = 1.3
        zs_1_margin = 0.001

        for (zs, indices) in zip(
                self.zs, self.sfit_files['63'].sfit_nob_indices):

            # sfit uses different calculations for z < 0.001 and z > 10.
            index_10 = (zs < 10.)
            index_001 = (zs < 0.001)
            # The sfit code is not accurate near 1.0.
            not_near_1 = (np.abs(zs - 1.) > 0.03)
            self.indices_mag_test.append(index_10 & ~index_001 & not_near_1)

            index_large = (zs > z_break)
            index_small = (zs <= z_break)
            self._indexes.append([index_large, index_small])

            # The sfit code is not accurate near 1.0.
            near_1 = (np.abs(zs - 1.) > zs_1_margin)
            self._indices_not_near_1.append(near_1)
            near_1_db = (zs < 0.88) | (zs > 1.1)
            self._indices_not_near_1_db.append(near_1_db)

    def test_get_magnification(self):
        for (nob_indices, mag_test_indices, gamma, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.gammas, self.mag_objs):
            mag = mag_obj.get_magnification()

            sfit_mag = self.sfit_files['61'].mag[nob_indices][mag_test_indices]
            b1_factor = (self.sfit_files['63'].amp[nob_indices][
                             mag_test_indices] *
                         self.sfit_files['63'].b1[nob_indices][
                             mag_test_indices] * gamma)
            sfit_mag = sfit_mag + b1_factor

            np.testing.assert_allclose(
                mag[mag_test_indices], sfit_mag, rtol=0.005)

    def test_magnification(self):
        for (nob_indices, mag_test_indices, gamma, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.gammas, self.mag_objs):

            with self.assertRaises(AttributeError):
                _ = mag_obj.magnification

            sfit_mag = self.sfit_files['61'].mag[nob_indices][mag_test_indices]
            b1_factor = (self.sfit_files['63'].amp[nob_indices][
                             mag_test_indices] *
                         self.sfit_files['63'].b1[nob_indices][
                             mag_test_indices] * gamma)
            sfit_mag += b1_factor

            mag_obj.get_magnification()
            np.testing.assert_allclose(
                mag_obj.magnification[mag_test_indices],
                sfit_mag, rtol=0.005)

    def _get_factor_b0(self, nob_indices):
        fspl_factor_b0 = (self.sfit_files['63'].amp[nob_indices] *
                          self.sfit_files['61'].db0[nob_indices])
        fspl_factor_b0 /= self.sfit_files['51'].a[3]  # rho
        fspl_factor_b0 += (self.sfit_files['62'].dAdu[nob_indices] *
                           self.sfit_files['63'].b0[nob_indices])

        return fspl_factor_b0

    def test_get_d_u_d_params(self):
        super().test_get_d_u_d_params()

    def test_get_d_A_d_u(self):
        self._get_d_A_d_u_1()
        self._get_d_A_d_u_2()

    def _get_d_A_d_u_1(self):
        """
        sfit returns: FSPL:
            61 dA/drho
            62 df/dparams, dAdu

        A_US = A_PS(u) * b0(z)
        dAdu_US = dA_PS(u) * b0(z) + A_PS(u) * db0(z)

        dAdu_FS = (damp*(b0-gamma*b1)+amp*(db0-gamma*db1)/rho)
         = dA*b0 - dA*gamma*b1 + A*db0/rho - A*gamma*db1/rho
         = dAdu_US - dA*gamma*b1 - A*gamma*db1/rho
        """
        rho = self.sfit_files['51'].a[3]

        for (nob_indices, mag_test_indices, mag_obj, gamma) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs,
                self.gammas):
            dA_du = mag_obj.get_d_A_d_u()[mag_test_indices]

            sfit_dA_du_US = (self.sfit_files['63'].b0[nob_indices][mag_test_indices] *
                             self.sfit_files['62'].dAdu[nob_indices][mag_test_indices])
            sfit_dA_du_US += (self.sfit_files['63'].amp[nob_indices][mag_test_indices] *
                              self.sfit_files['61'].db0[nob_indices][mag_test_indices] /
                              rho)
            np.testing.assert_allclose(dA_du, sfit_dA_du_US, rtol=0.015)

    def _get_d_A_d_u_2(self):
        """
        sfit returns: FSPL:
            61 dA/drho
            62 df/dparams, dAdu

        A_US = A_PS(u) * b0(z)
        dAdu_US = dA_PS(u) * b0(z) + A_PS(u) * db0(z)

        dAdu_FS = (damp*(b0-gamma*b1)+amp*(db0-gamma*db1)/rho)
         = dA*b0 - dA*gamma*b1 + A*db0/rho - A*gamma*db1/rho
         = dAdu_US - dA*gamma*b1 - A*gamma*db1/rho

        df_FS = fs * (dAdu_US - dA*gamma*b1 - A*gamma*db1/rho)
        """
        rho = self.sfit_files['51'].a[3]

        zip_ = zip(self.sfit_files['62'].sfit_nob_indices, self.sfit_files['51'].source_fluxes,
                   self.gammas, self.indices_mag_test, self.mag_objs)
        for (nob_indices, source_flux, gamma, mag_test_indices, mag_obj) in zip_:
            dA_du = mag_obj.get_d_A_d_u()[mag_test_indices]

            b1_term = (self.sfit_files['62'].dAdu * self.sfit_files['63'].b1)[nob_indices][mag_test_indices]
            b1_term += (self.sfit_files['63'].amp * self.sfit_files['61'].db1)[nob_indices][mag_test_indices] / rho
            b1_term *= gamma
            sfit_dA_du_US = self.sfit_files['62'].df[nob_indices][mag_test_indices] / source_flux + b1_term

            np.testing.assert_allclose(dA_du, sfit_dA_du_US, rtol=0.015)

    def test_get_d_A_d_params(self):
        """
        df/dparams = fs * dA/dparams (FSPL)

        dA/dparams (PSPL) = d_A_d_u * d_u_d_params[key]
        dA/dparams (FSPL) = d_u_d_params[key] * factor
            = factor * dA/dparams(PSPL) / dA_du

        dA/dparams (b0) = d_u_d_params * factor_b0
        dA/dparams (FSPL) = d_u_d_params * factor_b1

        dA/dparams (b0) = dA/dparams (FSPL) * factor_b1 / factor_b0
         = df/dparams * factor_b1 / factor_b0 / fs

        b0:
        factor = self.pspl_magnification * self.db0
        factor /= self.trajectory.parameters.rho
        factor += self.get_d_A_d_u() * self.b0

        b0, b1:
        factor = self.pspl_magnification * (self.db0 - self._gamma * self.db1)
        factor /= self.trajectory.parameters.rho
        factor += self.get_d_A_d_u() * (self.b0 - self._gamma * self.b1)
        """
        params = ['t_0', 'u_0', 't_E']

        for (nob_indices, source_flux, gamma, mag_test_indices, not_near_1,
             mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.gammas,
                self.indices_mag_test, self._indices_not_near_1_db,
                self.mag_objs):

            dA_dparam = mag_obj.get_d_A_d_params(params)
            fspl_factor_b1 = self._get_factor_b1(nob_indices, gamma)
            fspl_factor_b0 = self._get_factor_b0(nob_indices)

            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                sfit_df_dparam = self.sfit_files['62'].data[
                    'dfd{0}'.format(short_param)][nob_indices]

                sfit_dA_dparam = (sfit_df_dparam * fspl_factor_b0 /
                                  fspl_factor_b1 / source_flux)

                np.testing.assert_allclose(
                    dA_dparam[param][mag_test_indices & not_near_1],
                    sfit_dA_dparam[mag_test_indices & not_near_1], rtol=0.015)

    def test_get_d_A_d_rho(self):
        """
        d_A_d_rho = np.ones(len(self.trajectory.times))
        d_A_d_rho *= self.pspl_magnification
        d_A_d_rho *= -self.u_ / self.trajectory.parameters.rho**2
        d_A_d_rho *= (self.db0 - self._gamma * self.db1)

        dA_drho_b0 = db0 * dA_drho_b1 /(db0 - gamma*db1)
        """
        for (nob_indices, source_flux,  gamma, mag_test_indices,
             mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes, self.gammas,
                self.indices_mag_test, self.mag_objs):

            sfit_df_dparam = self.sfit_files['61'].data['dfdrho'][nob_indices]
            factor = self.sfit_files['61'].data['db0'][nob_indices]
            factor /= (self.sfit_files['61'].data['db0'][nob_indices] -
                       gamma * self.sfit_files['61'].data['db1'][nob_indices])
            sfit_dA_drho = factor * sfit_df_dparam / source_flux
            dAdrho = mag_obj.get_d_A_d_rho()

            np.testing.assert_allclose(
                dAdrho[mag_test_indices], sfit_dA_drho[mag_test_indices],
                rtol=0.015)

    def test_b0(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):

            np.testing.assert_allclose(
                mag_obj.b0[mag_test_indices],
                self.sfit_files['63'].b0[nob_indices][mag_test_indices],
                atol=0.0001)

    def test_db0(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):

            np.testing.assert_allclose(
                mag_obj.db0[mag_test_indices],
                self.sfit_files['61'].db0[nob_indices][mag_test_indices],
                atol=0.005)

    def test_z_(self):
        for (nob_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices, self.mag_objs):
            np.testing.assert_allclose(
                mag_obj.z_,
                self.sfit_files['63'].x[nob_indices] /
                self.sfit_files['51'].a[3],
                rtol=0.0001)


class TestFiniteSourceUniformGould94DirectMagnification(TestFiniteSourceUniformGould94Magnification):

    def setUp(self):
        TestFiniteSourceUniformGould94Magnification.setUp(self)

    def _set_mag_objs(self):
        for trajectory in self.trajectories:
            mag_obj = mm.FiniteSourceUniformGould94Magnification(
                trajectory=trajectory, direct=True)
            self.mag_objs.append(mag_obj)

    def test_get_d_A_d_u(self):
        # derivatives of B_0 not implemented for direct method.
        pass

    def test_db0(self):
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                _ = mag_obj.db0

    def test_get_d_A_d_params(self):
        """
        df/dparams = fs * dA/dparams (FSPL)
        """
        params = ['t_0', 'u_0', 't_E', 'rho']
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                mag_obj.get_d_A_d_params(params)

    def test_get_d_A_d_rho(self):
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                mag_obj.get_d_A_d_rho()


class TestFiniteSourceLDYoo04Magnification(TestFiniteSourceUniformGould94Magnification):

    def setUp(self):
        TestFiniteSourceUniformGould94Magnification.setUp(self)

    def _set_mag_objs(self):
        for (trajectory, gamma) in zip(self.trajectories, self.gammas):
            mag_obj = mm.FiniteSourceLDYoo04Magnification(
                trajectory=trajectory, gamma=gamma)
            self.mag_objs.append(mag_obj)

    def test_get_magnification(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):
            mag = mag_obj.get_magnification()

            np.testing.assert_allclose(
                mag[mag_test_indices],
                self.sfit_files['61'].mag[nob_indices][mag_test_indices],
                rtol=0.005)

    def test_magnification(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):

            with self.assertRaises(AttributeError):
                _ = mag_obj.magnification

            mag_obj.get_magnification()
            np.testing.assert_allclose(
                mag_obj.magnification[mag_test_indices],
                self.sfit_files['61'].mag[nob_indices][mag_test_indices],
                rtol=0.005)

    def test_get_d_u_d_params(self):
        super().test_get_d_u_d_params()

    def test_get_d_A_d_u(self):
        """
        Straight up reported in fort.62.
        """
        for (nob_indices, source_flux, mag_test_indices, mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.indices_mag_test, self.mag_objs):

            dA_du = mag_obj.get_d_A_d_u()
            sfit_dA_du = self.sfit_files['62'].df[nob_indices] / source_flux
            np.testing.assert_allclose(
                dA_du[mag_test_indices],
                sfit_dA_du[mag_test_indices], rtol=0.015)

    def test_get_d_A_d_params(self):
        """
        df/dparams = fs * dA/dparams (FSPL)
        """
        params = ['t_0', 'u_0', 't_E']

        for (nob_indices, source_flux, gamma, mag_test_indices, mag_obj) in zip(
                self.sfit_files['62'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.gammas, self.indices_mag_test, self.mag_objs):

            dA_dparam = mag_obj.get_d_A_d_params(params)

            for j, param in enumerate(params):
                short_param = param.replace('_', '')
                sfit_df_dparam = self.sfit_files['62'].data[
                    'dfd{0}'.format(short_param)][nob_indices]
                sfit_dA_dparam = sfit_df_dparam / source_flux
                np.testing.assert_allclose(
                    dA_dparam[param][mag_test_indices],
                    sfit_dA_dparam[mag_test_indices], rtol=0.015)

    def test_get_d_A_d_rho(self):
        for (nob_indices, source_flux, mag_test_indices, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.sfit_files['51'].source_fluxes,
                self.indices_mag_test, self.mag_objs):

            sfit_df_dparam = self.sfit_files['61'].data['dfdrho'][nob_indices]
            sfit_dA_dparam = sfit_df_dparam / source_flux
            dAdrho = mag_obj.get_d_A_d_rho()

            np.testing.assert_allclose(
                dAdrho[mag_test_indices], sfit_dA_dparam[mag_test_indices],
                rtol=0.015)

    def test_b1(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['63'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):

            np.testing.assert_allclose(
                mag_obj.b1[mag_test_indices],
                self.sfit_files['63'].b1[nob_indices][mag_test_indices],
                atol=0.0001)

    def test_db1(self):
        for (nob_indices, mag_test_indices, mag_obj) in zip(
                self.sfit_files['61'].sfit_nob_indices,
                self.indices_mag_test, self.mag_objs):

            np.testing.assert_allclose(
                mag_obj.db1[mag_test_indices],
                self.sfit_files['61'].db1[nob_indices][mag_test_indices],
                atol=0.001)

    def test_gamma(self):
        for (gamma, mag_obj) in zip(self.gammas, self.mag_objs):
            np.testing.assert_almost_equal(gamma, mag_obj.gamma)


class TestFiniteSourceLDYoo04DirectMagnification(TestFiniteSourceLDYoo04Magnification):

    def setUp(self):
        TestFiniteSourceLDYoo04Magnification.setUp(self)

    def _set_trajectories(self):
        # Take only N_MAX indices from within t_0 +- 3*t_star to reduce runtime
        n_max = 10

        index_1 = self.sfit_files['63'].t > (
                self.parameters.t_0 - 3. * self.parameters.t_star)
        index_2 = self.sfit_files['63'].t < (
                self.parameters.t_0 + 3. * self.parameters.t_star)
        index = index_1 & index_2

        for i, nob_indices in enumerate(self.sfit_files['63'].sfit_nob_indices):
            new_nob_indices = nob_indices & index
            if np.sum(new_nob_indices) > n_max:
                n_skip = np.floor(np.sum(new_nob_indices) / n_max).astype(int)
                mask = np.zeros(len(nob_indices), dtype=bool)
                mask[::n_skip] = True
                new_nob_indices = new_nob_indices & mask

            for f_sixty in ['61', '62', '63']:
                self.sfit_files[f_sixty].sfit_nob_indices[i] = new_nob_indices

        TestFiniteSourceLDYoo04Magnification._set_trajectories(self)

    def _set_mag_objs(self):
        for (trajectory, gamma) in zip(self.trajectories, self.gammas):
            mag_obj = mm.FiniteSourceLDYoo04Magnification(
                trajectory=trajectory, gamma=gamma, direct=True)
            self.mag_objs.append(mag_obj)

    def test_get_d_A_d_u(self):
        # derivatives of B_0 not implemented for direct method.
        pass

    def test_db0(self):
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                _ = mag_obj.db0

    def test_db1(self):
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                _ = mag_obj.db1

    def test_get_d_A_d_params(self):
        """
        df/dparams = fs * dA/dparams (FSPL)
        """
        params = ['t_0', 'u_0', 't_E', 'rho']
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                mag_obj.get_d_A_d_params(params)

    def test_get_d_A_d_rho(self):
        for mag_obj in self.mag_objs:
            with self.assertRaises(NotImplementedError):
                mag_obj.get_d_A_d_rho()
