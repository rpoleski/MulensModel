# import matplotlib
# matplotlib.use('Agg')
import sys
import yaml
import numpy as np
from astropy import units as u

from ulens_model_fit import UlensModelFit

from modelIMF import ModelIMF


kappa = 8.144  # [mas / M_Sun]
const_4_74 = 4.7405  # v_t [km/s] = 4.74 * mu [mas/yr] / pi [mas]


class UlensModelFitWithGalacticModel(UlensModelFit):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.model_IMF = ModelIMF()
        self.user_constants = None

    def _set_default_parameters(self):
        """
        Extend the set of available parameters
        """
        super()._set_default_parameters()
        self._other_parameters = ['M_l', 'D_l', 'mu_s_N', 'mu_s_E', 'mu_rel']
        self._latex_conversion_other = dict(
            M_l='M_{\\rm lens}',
            D_l='D_{\\rm lens}',
            mu_s_N='\\mu_{{\\rm s}_N}',
            mu_s_E='\\mu_{{\\rm s}_E}',
            mu_rel="\\mu_{\\rm rel}")

    def _get_samples_for_triangle_plot(self):
        """
        Prepare a cube of samples for triangle plot
        """
        parameters_to_add = "M_l v_rel pi_rel".split()

        n_add = len(parameters_to_add)
        out = np.zeros(np.array(self._samples_flat.shape) + (0, n_add))
        out[:, :-n_add] = self._samples_flat

        # TODO: extract the lines below from existing info
        indexes = [2, 3, 4, 5, 6, 7]
        keys = "t_E pi_E_N pi_E_E mu_s_N mu_s_E mu_rel".split()
        kwargs = {k: self._samples_flat[:, i] for (k, i) in zip(keys, indexes)}

        kwargs['D_s'] = self.D_s
        params = self._get_parameters(**kwargs)

        for (i, parameter) in enumerate(parameters_to_add):
            out[:, -n_add+i] = params[parameter]

        # TODO - this should be in settings:
        # here we change mass -> log10(MASS) for triangle plot
        out[:, -3] = np.log10(out[:, -3])

        return out

    def _get_labels_for_triangle_plot(self):
        """
        provide list of labels to be used by triangle plot
        """
        # TODO - this should be in settings:
        out = self._fit_parameters_latex + [
            "$\\log_{10}(M_l)$", "$v_{\\rm rel}$", "$\\pi_{\\rm rel}$"]
        return out

# BELOW IS GALACTIC MODEL
    def set_user_constants(self):
        """
        Set coordinates-depending constants

        All velocities are in km/s
        """
        self.V_Sun = 232.2
        self.W_Sun = 7.3

        # These velocities are taken from Batista et al. (2011)
        self.v_y_disk = 220.
        self.v_z_disk = 0.
        self.sigma_y_disk = 30.
        self.sigma_z_disk = 20.
        self.sigma_y_bulge = 100.
        self.sigma_z_bulge = 100.

        self.R_Sun = 8.178  # Based on:
        # https://ui.adsabs.harvard.edu/abs/2019A&A...625L..10G

        self.D_s = 8.5  # assumed source distance
        self.bulge_angle = 20.0 * np.pi / 180.0  # assumed bulge angle is 20deg

        ra_G_rad = 192.85948 * np.pi / 180.
        dec_G_rad = 27.12825 * np.pi / 180.

        l_rad = self._model.coords.galactic.l.to(u.rad).value
        b_rad = self._model.coords.galactic.b.to(u.rad).value
        ra_rad = self._model.coords.ra.to(u.rad).value
        dec_rad = self._model.coords.dec.to(u.rad).value

        self.const_x = np.cos(b_rad) * np.sin(l_rad)
        self.const_y = np.cos(b_rad) * np.cos(l_rad)
        self.const_z = np.sin(b_rad)

        C_1 = (np.sin(dec_G_rad) * np.cos(dec_rad) -
               np.cos(dec_G_rad) * np.sin(dec_rad) * np.cos(ra_rad - ra_G_rad))
        self.C_1 = C_1 / np.cos(b_rad)
        self.C_2 = np.cos(dec_G_rad) * np.sin(ra_rad-ra_G_rad) / np.cos(b_rad)

        t_0_par = self._model.parameters.t_0_par
        v_Earth_perp = self._model.coords.v_Earth_projected(t_0_par)
        self.v_Earth_perp_N = v_Earth_perp[0]
        self.v_Earth_perp_E = v_Earth_perp[1]

        self.user_constants = "DONE"

    def _get_ln_probability_for_other_parameters(self):
        """
        Calculate log(event probability) as a function of added parameters.
        Here we use event rate for that purpose.
        """
        if self.user_constants is None:
            self.set_user_constants()

        t_E = self._model.parameters.t_E
        pi_E_N = self._model.parameters.pi_E_N
        pi_E_E = self._model.parameters.pi_E_E
        mu_s_N = self._other_parameters_dict['mu_s_N']
        mu_s_E = self._other_parameters_dict['mu_s_E']
        mu_rel = self._other_parameters_dict['mu_rel']

        params = self._get_parameters(
            t_E=t_E, pi_E_N=pi_E_N, pi_E_E=pi_E_E,
            mu_s_N=mu_s_N, mu_s_E=mu_s_E, mu_rel=mu_rel, D_s=self.D_s)

        ln_mass_probability = self.model_IMF.get_ln_relative_probability(
            params['M_l'])
        if ln_mass_probability == -np.inf:
            return -np.inf

        rho_disk = self.rho_disk(x=params['x'], y=params['y'], z=params['z'])
        rho_bulge = self.rho_bulge(x=params['x'], y=params['y'],
                                   z=params['z'])

        rho_disk_source = self.rho_disk(x=params['x_s'], y=params['y_s'],
                                        z=params['z_s'])
        rho_bulge_source = self.rho_bulge(x=params['x_s'], y=params['y_s'],
                                          z=params['z_s'])

        difference_y = (params['V_l_l'] - self.v_y_disk) / self.sigma_y_disk
        difference_z = (params['V_l_b'] - self.v_z_disk) / self.sigma_z_disk
        probability_disk = np.exp(-0.5*(difference_y**2 + difference_z**2))
        probability_disk /= self.sigma_y_disk * self.sigma_z_disk

        temp_1 = (params['V_l_l'] / self.sigma_y_bulge)**2
        temp_2 = (params['V_l_b'] / self.sigma_z_bulge)**2
        temp_3 = self.sigma_y_bulge * self.sigma_z_bulge
        probability_bulge = np.exp(-0.5*(temp_1 + temp_2)) / temp_3

        # source:
        diff_y = (params['V_s_l'] - self.v_y_disk) / self.sigma_y_disk
        diff_z = (params['V_s_b'] - self.v_z_disk) / self.sigma_z_disk
        probability_disk_source = np.exp(-0.5*(diff_y**2 + diff_z**2))
        probability_disk_source /= self.sigma_y_disk * self.sigma_z_disk

        temp_1 = (params['V_s_l'] / self.sigma_y_bulge)**2
        temp_2 = (params['V_s_b'] / self.sigma_z_bulge)**2
        temp_3 = self.sigma_y_bulge * self.sigma_z_bulge
        probability_bulge_source = np.exp(-0.5*(temp_1 + temp_2)) / temp_3

        ln_rate = ln_mass_probability + np.log(params['M_l'])
        ln_rate += 4. * np.log(params['D_l'] * mu_rel)
        ln_rate += np.log(t_E)
        ln_rate -= np.log(params['pi_E'])
        ln_rate += np.log(params['D_l']**2 * (
            rho_disk * probability_disk + rho_bulge * probability_bulge))

        ln_rate += np.log(rho_disk_source * probability_disk_source +
                          rho_bulge_source * probability_bulge_source)

        return ln_rate

    def rotate_xy_to_bulge_axis(self, x, y):
        """
        rotate coordinates so that they are along bulge axis
        """
        x_ = x*np.sin(self.bulge_angle) + y*np.cos(self.bulge_angle)
        y_ = x*np.cos(self.bulge_angle) - y*np.sin(self.bulge_angle)
        return (x_, y_)

    def rho_bulge(self, x, y, z):
        """
        bulge density based on
        Dwek 1995 1995ApJ...445..716D

        Bulge cut-off radius is taken from Robin (2003) ???
        """
        (xp, yp) = self.rotate_xy_to_bulge_axis(x, y)
        r = np.sqrt(x**2 + y**2)
        rs_4 = ((xp/1.58)**2 + (yp/0.62)**2)**2 + (z/0.43)**4
        rho = 1.23 * np.exp(-0.5*np.sqrt(rs_4))
        if r > 2.4:
            rho *= np.exp(-0.5*((r-2.4)/0.5)**2)
        return rho

    def rho_disk(self, x, y, z):
        """
        Disk density based on Batista+11 Table 2
        """
        beta = 0.381
        R = np.sqrt(x**2 + y**2)
        rho = (1.0-beta)*np.exp(-abs(z)/0.156) + beta*np.exp(-abs(z)/0.439)
        rho = rho * 1.07 * np.exp(-R/2.75)
        if R < 1.0:
            rho = 0.0
        return rho

    def _get_parameters(self, t_E, pi_E_N, pi_E_E,
                        mu_s_N, mu_s_E, mu_rel, D_s):
        """
        Calculate physical model parameters
        """
        t_E_year = t_E / 365.25

        out = dict()
        out['pi_E'] = np.sqrt(pi_E_N**2 + pi_E_E**2)
        out['theta_E'] = mu_rel * t_E_year
        out['pi_rel'] = out['theta_E'] * out['pi_E']
        out['M_l'] = out['theta_E'] / (kappa * out['pi_E'])
        out['D_l'] = 1. / (out['pi_rel'] + 1. / D_s)

        scale = mu_rel / out['pi_E']
        out['mu_l_N'] = mu_s_N + scale * pi_E_N
        out['mu_l_E'] = mu_s_E + scale * pi_E_E
        mu_l_NE_helio = self.shift_velocity_from_geo_to_helio(
            out['mu_l_N'], out['mu_l_E'], out['pi_rel'])
        out['mu_l_N_helio'] = mu_l_NE_helio[0]
        out['mu_l_E_helio'] = mu_l_NE_helio[1]
        mu_l_lb_helio = self.rotate_pm_radec_to_lb(out['mu_l_E_helio'],
                                                   out['mu_l_N_helio'])
        out['mu_l_l_helio'] = mu_l_lb_helio[0]
        out['mu_l_b_helio'] = mu_l_lb_helio[1]
        (out['V_l_l'], out['V_l_b']) = self.get_velocity_from_pm(
            out['mu_l_l_helio'], out['mu_l_b_helio'], out['D_l'])

        mu_s_NE_helio = self.shift_velocity_from_geo_to_helio(mu_s_N, mu_s_E,
                                                              out['pi_rel'])
        out['mu_s_N_helio'] = mu_s_NE_helio[0]
        out['mu_s_E_helio'] = mu_s_NE_helio[1]
        mu_s_lb_helio = self.rotate_pm_radec_to_lb(out['mu_s_E_helio'],
                                                   out['mu_s_N_helio'])
        out['mu_s_l_helio'] = mu_s_lb_helio[0]
        out['mu_s_b_helio'] = mu_s_lb_helio[1]
        (out['V_s_l'], out['V_s_b']) = self.get_velocity_from_pm(
            out['mu_s_l_helio'], out['mu_s_b_helio'], D_s)

        out['x'] = out['D_l'] * self.const_x
        out['y'] = self.R_Sun - out['D_l'] * self.const_y
        out['z'] = out['D_l'] * self.const_z

        out['x_s'] = D_s * self.const_x
        out['y_s'] = self.R_Sun - D_s * self.const_y
        out['z_s'] = D_s * self.const_z

        out['v_rel'] = const_4_74 * mu_rel / out['pi_rel']
        return out

    def shift_velocity_from_geo_to_helio(self, velocity_N, velocity_E, pi_rel):
        """
        Calculate heliocentric proper motion based on geocentric
        """
        vel_N_helio = velocity_N + self.v_Earth_perp_N * pi_rel / const_4_74
        vel_E_helio = velocity_E + self.v_Earth_perp_E * pi_rel / const_4_74
        return (vel_N_helio, vel_E_helio)

    def rotate_pm_radec_to_lb(self, mu_E, mu_N):
        """
        rotate proper motion from RA/Dec to galactic l/b coords
        """
        mu_l = self.C_1 * mu_E + self.C_2 * mu_N
        mu_b = -self.C_2 * mu_E + self.C_1 * mu_N
        return (mu_l, mu_b)

    def get_velocity_from_pm(self, mu_l, mu_b, distance):
        """
        Calculate velocity based on proper motion and distance
        """
        V_l = const_4_74 * mu_l * distance + self.V_Sun
        V_b = const_4_74 * mu_b * distance + self.W_Sun
        return (V_l, V_b)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = UlensModelFitWithGalacticModel(**settings)

    ulens_model_fit.run_fit()
