"""Code for fitting microlensing models with MIST evolutionary tracks using implemented in Brutus pagege.
for more details see: Mróz, M.J. el al. 2025, https://doi.org/10.1051/0004-6361/202660534"""
import os
import sys
import yaml
import numpy as np
import math
import shutil
from copy import copy
import matplotlib.pyplot as plt
try:
    from brutus import EEPTracks
    from brutus.core import StarEvolTrack
except ImportError:
    print("Warning: Brutus is not installed. Please install Brutus to use evolutionary tracks fitting.")
    print("see: https://github.com/joshspeagle/brutus")

import MulensModel as mm
from ulens_model_fit import UlensModelFit, import_failed

R_sun_over_au = 4.650467260962158  # [1000 * R_sun/au]


class UlensModelFitEvolutionaryTracks(UlensModelFit):
    def __init__(self, photometry_files, evolutionary_tracks=None, starting_parameters=None, prior_limits=None,
                 model=None, fixed_parameters=None, extra_parameters=None, min_values=None, max_values=None,
                 fitting_parameters=None, fit_constraints=None, plots=None, other_output=None, fit_method=None):
        self._plot_track_file = plots.pop('track', {}).get('file') if plots is not None else None
        self._parce_brutus_settings(evolutionary_tracks)
        super().__init__(photometry_files=photometry_files, starting_parameters=starting_parameters,
                         prior_limits=prior_limits, model=model, fixed_parameters=fixed_parameters,
                         extra_parameters=extra_parameters, min_values=min_values, max_values=max_values,
                         fitting_parameters=fitting_parameters, fit_constraints=fit_constraints,
                         plots=plots, other_output=other_output, fit_method=fit_method)

    def _parce_brutus_settings(self, evolutionary_tracks):
        """Read and check the Brutus settings from the YAML file."""
        self._brutus_data = None
        self._photometric_map_file = None
        self._loga_max = 13.  # Kpc
        self._eep_max = 808
        self._brutus_tol = 1e-5
        self._sigma_theta_E_bound = 0.1  # mas
        if evolutionary_tracks is not None:
            self._brutus_data = evolutionary_tracks.get('brutus_data', None)
            self._loga_max = float(evolutionary_tracks.get('loga_max', 13.))
            self._eep_max = float(evolutionary_tracks.get('eep_max', 808))
            self._brutus_tol = float(evolutionary_tracks.get('tol', 1e-5))
            self._sigma_theta_E_bound = float(evolutionary_tracks.get('sigma_theta_E_bound', 0.1))
            self._photometric_map_file = evolutionary_tracks.get('photometric_map_file', None)

        if self._brutus_data is None:
            self._brutus_data = os.environ.get('BRUTUS_DATA', None)
        if self._brutus_data is None:
            raise ValueError("BRUTUS_DATA environment variable is not set. Please set it to the path" +
                             " of the Brutus grids directory or define it in the input YAML file.")

        if self._photometric_map_file is None:
            print(
                'Warning: photometric_map_file is not set. The color-magnitude '
                'diagram (CMD) will include only the evolutionary tracks and will '
                'not show field stars.'
            )
            print(
                'To include field stars in the CMD, set photometric_map_file to a '
                'file containing field-star magnitudes in the same photometric '
                'bands used by the MIST models.'
            )

    def imf_lnprior(self, mass, alpha_low=1.16, alpha_high=2.32, mass_break=0.9):
        """
        Apply a Kroupa-like broken IMF prior over the provided initial mass grid.
        parameters base on
        Koshimoto https://iopscience.iop.org/article/10.3847/1538-4357/ac07a8/pdf
        Parameters
        ----------
        mass : `float`
             initial mass (solar units) the IMF will be evaluated over.

        alpha_low : float, optional
            Power-law slope for the low-mass component of the IMF.
            Default is `1.16`.

        alpha_high : float, optional
            Power-law slope for the high-mass component of the IMF.
            Default is `2.32`.

        mass_break : float, optional
            The mass where we transition from `alpha_low` to `alpha_high`.
            Default is `0.9`.
        """

        # Initialize log-prior.
        lnprior = -np.inf

        # Low mass.
        if (mass <= mass_break) & (mass > 0.08):
            lnprior = -alpha_low * np.log(mass)

        # High mass.
        if mass > mass_break:
            lnprior = (-alpha_high * np.log(mass)+(alpha_high - alpha_low) * np.log(mass_break))

        # Compute normalization.
        norm_low = mass_break**(1.-alpha_low)/(alpha_high - 1.)
        norm_high = 0.08**(1.-alpha_low)/(alpha_low-1.)  # H-burning limit
        norm_high -= mass_break**(1.-alpha_low)/(alpha_low-1.)
        norm = norm_low + norm_high

        return lnprior - np.log(norm)

    def feh_lnprior(self, feh, feh_mean=-0.2, feh_sigma=0.3):
        """
        Log-prior for the metallicity in a given component of the galaxy.

        Parameters
        ----------
        feh : `~numpy.ndarray` of shape (N)
            The metallicities of the corresponding models.

        feh_mean : float, optional
            The mean metallicity. Default is `-0.2`.

        feh_sigma : float, optional
            The standard deviation in the metallicity. Default is `0.3`.

        Returns
        -------
        logp : `~numpy.ndarray` of shape (N)
            The corresponding normalized ln(probability).
        """

        # Compute log-probability.
        chi2 = (feh_mean - feh)**2 / feh_sigma**2  # chi2
        lnorm = np.log(2. * np.pi * feh_sigma**2)  # normalization
        lnprior = -0.5 * (chi2 + lnorm)
        return lnprior

    def _get_theta_star_from_isochrone(self, track):
        """
        Get the theta_star from the isochrones models.
        """
        D_S = self._other_parameters_dict['D_S']
        radius = 10**track['logr']
        return R_sun_over_au * radius / D_S  # [1000 * au ]

    def _get_theta_E_from_rho_and_isochrone(self):
        """
        Calculate the Einstein radius from the isochrone model and rho
        """
        # only primary source
        track = self._source_parameters[0]
        rho = self._model.parameters.parameters.get('rho_1', self._model.parameters.parameters.get('rho', 0.))
        # D_S = self._other_parameters_dict['D_S']
        self._theta_E_rho = self._get_theta_star_from_isochrone(track) / rho   # [mas]
        # print(f'theta_E_rho={theta_E_rho}, rho={rho}, D_S={D_S}')
        return self._theta_E_rho

    def _get_theta_E_from_xallarap_and_isochrone(self):
        """
        Calculate the Einstein radius from the isochrone model and xallarap orbit
        """
        a_physical = self._get_a_physical()  # [au]
        D_S = self._other_parameters_dict['D_S']
        # a_einstein = self._model.parameters.parameters['xi_semimajor_axis']  # [1/theta_E]
        theta_E_xi = 1000 * a_physical / \
            self._model.parameters.parameters['xi_semimajor_axis'] / D_S  # [1000 *  au / pc = mas]
        # print(f'a_physical={a_physical}, D_S={D_S}, a_einstein={a_einstein}, theta_E_xi={theta_E_xi}')
        return theta_E_xi

    def _get_a_physical(self):
        """Get semi-major axis in physical units [au] from xallarap orbit and isochrone model.
        """

        # M_S_1 = self._other_parameters_dict['M_S']
        M_S_1 = self._source_parameters[0]['mass']
        if 'q_source' in self._fit_parameters:
            q_S = self._model.parameters.parameters['q_source']
        else:
            q_S = self._other_parameters_dict['q_S']
        xi_period = self._model.parameters.parameters['xi_period'] / 365.25
        return np.cbrt((M_S_1 * (1. + q_S)) / (xi_period**2.))  # [au]

    def _set_default_parameters(self):
        """
        Extend the set of available parameters
        """
        self._reparametrized_MM_parameters_values = None
        self._reparametrized_MM_parameters = None
        super()._set_default_parameters()
        self._other_parameters += [
            'M_ini_S', 'M_ini_S_1', 'M_ini_S_2', 'EEP_S', 'EEP_S_1', 'EEP_S_2', 'age_S', 'feh_S', 'D_S', 'AV_S', 'q_S']
        self._latex_conversion_other.update({
            'M_ini_S': 'M_{\\rm{ini, S}}',
            'M_ini_S_1': 'M_{\\rm{ini, S, 1}}',
            'M_ini_S_2': 'M_{\\rm{ini, S, 2}}',
            'EEP_S': 'EEP_{\\rm{S}}',
            'EEP_S_1': 'EEP_{\\rm{S, 1}}',
            'EEP_S_2': 'EEP_{\\rm{S, 2}}',
            'age_S': 'age_{\\rm{S}}',
            'feh_S': 'feh_{\\rm{S}}',
            'D_S': 'D_{\\rm{S}}',
            'AV_S': 'A_{\\rm{V,S}}',
            'q_S': 'q_{\\rm{S}}',
        })
        self._parameters_star_aux = ['mass', 'logl', 'loga', 'logt', 'logr', 'logg', 'feh_surf', 'afe_surf']

    def _setup_brutus(self, bandpass=None):
        """
        Set up the Brutus evolutionary track interpolation.
        """
        self._other_parameters_prior = None
        if bandpass is not None:
            self._MIST_bandpass = bandpass
        else:
            self._set_bandpass()
        print(f'Using bands: {self._MIST_bandpass}')
        brutus_data = self._brutus_data
        gridfile = os.path.join(brutus_data,  'grid_mist_v9.h5')
        mistfile = os.path.join(brutus_data, 'MIST_1.2_EEPtrk.h5')
        nnfile = os.path.join(brutus_data, 'nn_c3k.h5')

        self._tracks = EEPTracks(mistfile=mistfile, predictions=[
            "mass", "logl", "loga", "logt", "logr", "logg", "feh_surf", "afe_surf"])
        self._star_track = StarEvolTrack(tracks=self._tracks, filters=self._MIST_bandpass, nnfile=nnfile)
        print(
            f'Brutus setup complete. Using grid file: {gridfile}, MIST file: {mistfile}, NN file: {nnfile},' +
            f' filters: {self._MIST_bandpass}')

        self._set_q_source_priors()
        self._set_other_parameters_priors()
        if 'q_source' in self._fit_parameters_unsorted:
            # trick to avoid error of MulensModel check, In the fitting process it will be set from ET models
            self._fixed_parameters['rho_2'] = 0.1

    def _set_bandpass(self):
        """
        Set the bandpass for the photometry files, witch is used for the Mist models.
        If not set dataset would not be used for evolutionary track fit.
        """
        self._MIST_bandpass = []
        for file in self._photometry_files:
            print(f'Parsing bandpass for file: {file}')
            MIST_band = file.pop('MIST_bandpass', None)
            if MIST_band is not None:
                if MIST_band in ['I', 'V']:
                    MIST_band = 'Bessell_' + MIST_band
                self._MIST_bandpass.append(MIST_band)

    def _set_q_source_priors(self):
        """
        Set the prior for the source mass ratio.
        """
        self._q_S = True
        self._q_source_prior = False
        self._q_S = False
        self._q_S_model_prior = False  # when q_source know from normal models defined by fix value in yaml file
        if 'rho_1' in self._fit_parameters_unsorted and 'q_source' in self._fit_parameters_unsorted:
            self._use_theta_E_bound = True
        else:
            self._use_theta_E_bound = False

    def _set_other_parameters_priors(self):
        """
        Set the priors for the other parameters.
        """
        # optional priors for other parameters
        self._other_parameters_prior = {
            # 'M_int_S':  self.imf_lnprior,
            # 'feh_S': self.feh_lnprior,
        }
        if self._q_S_model_prior:
            # Gaussian prior for q_S
            mu = copy(self._fixed_parameters['q_source'])
            sigma = 0.001

            def _q_gauss_prior(value, mu=mu, sigma=sigma):
                diff = value - mu
                return -0.5*(diff/sigma)**2 - math.log(math.sqrt(2*np.pi)*sigma)
            self._other_parameters_prior['q_S'] = _q_gauss_prior

    def _get_ln_probability_for_other_parameters(self):
        """
        We have to define this function and it has to return a float but in
        this case the change of ln_prob is coded in _ln_like().
        """
        out = 0.
        # out = super()._get_ln_probability_for_other_parameters()

        if self._other_parameters_prior is not None:
            for key, prior in self._other_parameters_prior.items():
                if key in self._other_parameters_dict:
                    value = self._other_parameters_dict[key]
                    out += prior(value)
                    # print('other prior', key, value, out)
        if self._use_theta_E_bound and 'q_source' in self._fit_parameters:
            out += self._theta_E_bound()
        return out

    def _theta_E_bound(self):
        """
        Check if theta_E from xallarap and from rho is consistent.
        """
        sigma = self._sigma_theta_E_bound
        out = 0.
        theta_E_rho = self.sources_dict['source_1']['theta_E_rho']
        theta_E_xallarap = self.sources_dict['source_2']['theta_E_xi']
        # print(f'theta_E_rho={theta_E_rho}, theta_E_xallarap={theta_E_xallarap}')
        out += self._get_ln_prior_for_1_parameter(theta_E_rho, ['gauss', theta_E_xallarap, sigma])
        # print(f'theta_E_rho={theta_E_rho}, theta_E_xallarap={theta_E_xallarap}, out={out}')
        return out

    def _get_seds(self, theta, parameters_star, EEP_S=None, loga_max=None):
        """when q_source is fitted it predicts the secondary EEP """
        parameters = dict(zip(self._fit_parameters, theta))
        smf = parameters.get('q_source', 0.0)

        mass_ini_S = parameters_star['M_ini_S']
        feh_S = parameters_star['feh_S']
        distance_S = parameters_star['D_S']
        AV_S = parameters_star['AV_S']
        if EEP_S is None:
            EEP_S = parameters_star['EEP_S']
        if loga_max is None:
            loga_max = self._loga_max
        mags, source_parameters1, source_parameters2 = self._star_track.get_seds(
            mini=mass_ini_S, feh=feh_S, eep=EEP_S,  av=AV_S, dist=distance_S,
            smf=smf, loga_max=loga_max, eep_binary_max=self._eep_max, tol=self._brutus_tol, combine_seds=False)
        # print(f'mags={mags}, source_parameters1={source_parameters1}, source_parameters2={source_parameters2}')
        mags, failed = self._parce_brutus_mags(mags)

        if np.all([np.isnan(mags[key]).all() for key in mags]):
            failed = True
            # print(mags, 'failed all nan')
        if self._model.n_sources == 2:
            return (mags, [source_parameters1, source_parameters2], failed)
        else:
            return (mags, [source_parameters1], failed)

    def _parce_brutus_mags(self, mags):
        """
        Check if the Brutus magnitudes are valid.
        """
        # print(f'Brutus magnitudes before check: {mags}')
        mags_dict = {}
        failed = False
        try:
            if self._model.n_sources == 1:
                if mags is np.nan:
                    mags = [np.nan for i in range(len(self._MIST_bandpass))]

                if len(mags) != len(self._MIST_bandpass):
                    raise ValueError(f'Brutus returned {len(mags)} magnitudes,' +
                                     f' but {len(self._MIST_bandpass)} were expected.')
                mags = [[m] for m in mags]
                mags_dict = dict(zip(self._MIST_bandpass, mags))
            if self._model.n_sources == 2:
                mags_dict = {key: [mags[j][i] for j in range(
                    self._model.n_sources)] for i, key in enumerate(self._MIST_bandpass)}
        except Exception as e:
            print(f'Error while parsing Brutus magnitudes: {e}')
            failed = True
        # print(f'Brutus magnitudes: {mags_dict}')
        return mags_dict, failed

    def _ln_like(self, theta):
        """
        likelihood function
        """
        self._set_model_parameters(theta)
        parameters_star = self._other_parameters_dict

        # print(parameters_star, '\n')
        (mags, self._source_parameters, failed) = self._get_seds(theta, parameters_star)
        # print(mags)
        if failed:
            return -np.inf
        self._set_q_source()
        self._set_sources_parameters()
        if 'q_source' in self._fit_parameters and 'rho_2' not in self._fit_parameters:
            self._set_rho_2_from_isochrone()

        fix_source_flux = self._get_fix_source_flux(mags)
        # print(f'fix_source_flux: {fix_source_flux}\n')
        fix_blend_flux = self._check_blend_fluxes(fix_source_flux)
        # print(f'fix_blend_flux: {fix_blend_flux}\n')
        self._event = mm.Event(self._datasets, self._model, fix_source_flux=fix_source_flux,
                               fix_blend_flux=fix_blend_flux)
        # print(f'self._event={self._event}\n')
        self._event.sum_function = 'numpy.sum'
        # self._set_n_fluxes()

        chi2 = self._event.get_chi2()
        out = -0.5 * chi2

        if self._print_model:
            self._print_current_model(theta, chi2)

        if self._task == 'fit' and len(self._other_parameters_dict) > 0:
            # print('other')
            out += self._get_ln_probability_for_other_parameters()

        if self._task == 'fit' and len(self._errorbars_parameters_dict) > 0:
            out += self._ln_prob_errors()
        return out

    def _set_sources_parameters(self):
        """
        Set source parameters based on isochrone models.
        """
        self.sources_dict = {}
        self.sources_values = []
        for i, track in enumerate(self._source_parameters):
            self.sources_dict[f'source_{i+1}'] = {}
            self.sources_dict[f'source_{i+1}']['theta_star'] = self._get_theta_star_from_isochrone(track)
            self.sources_values.append(self.sources_dict[f'source_{i+1}']['theta_star'])
            if i == 0:
                self.sources_dict[f'source_{i+1}']['theta_E_rho'] = self._get_theta_E_from_rho_and_isochrone()  # [mas]
                self.sources_values.append(self.sources_dict[f'source_{i+1}']['theta_E_rho'])
            if i == 1 and 'q_source' in self._fit_parameters:
                self.sources_dict[f'source_{i+1}']['theta_E_xi'] = self._get_theta_E_from_xallarap_and_isochrone()
                self.sources_values.append(self.sources_dict[f'source_{i+1}']['theta_E_xi'])
            for key in self._parameters_star_aux:
                self.sources_dict[f'source_{i+1}'][key] = track[key]
                self.sources_values.append(track[key])

    def _set_rho_2_from_isochrone(self):
        """
        Set rho_2 based on physical radius for secondary source and
        theta_E calculated from primary source physical radius and its rho.
        """
        theta_star_2 = self.sources_dict['source_2']['theta_star']
        theta_E_rho_1 = self.sources_dict['source_1']['theta_E_rho']
        self.rho_2 = theta_star_2 / theta_E_rho_1
        setattr(self._model.parameters, 'rho_2', self.rho_2)

    def _set_q_source(self):
        if 'q_source' in self._fit_parameters:
            q = self._model.parameters.parameters['q_source']
            self._other_parameters_dict['q_S'] = q
        if self._model.n_sources == 2:
            q = self._source_parameters[1]['mass'] / self._source_parameters[0]['mass']
            self._fixed_parameters['q_source'] = q
            self._other_parameters_dict['q_S'] = q
            setattr(self._model.parameters, 'q_source', q)

    def _get_fix_source_flux(self, mags):
        fix_source_flux = {}
        if self._model.n_sources == 1:
            for i, dataset in enumerate(self._datasets):
                if self._ET_fit[i]:
                    fix_source_flux[dataset] = self._get_flux_from_mag_safe(
                        mags[self._MIST_bandpass[i]][0])
        if self._model.n_sources == 2:
            for i, dataset in enumerate(self._datasets):
                if self._ET_fit[i]:
                    fix_source_flux[dataset] = [
                        self._get_flux_from_mag_safe(mags[self._MIST_bandpass[i]][0]),
                        self._get_flux_from_mag_safe(mags[self._MIST_bandpass[i]][1])]
                    # print('fix_source_flux ', dataset.plot_properties['label'], fix_source_flux[dataset])
        return fix_source_flux

    def _check_blend_fluxes(self, fix_source_flux):
        fix_blend_flux = {}
        # if self._model.n_sources == 1:
        #     for i, dataset in enumerate(self._datasets):
        #         if fix_source_flux[dataset] == 0.:
        #             fix_blend_flux[dataset] = 0.
        # if self._model.n_sources == 2:
        #     for i, dataset in enumerate(self._datasets):
        #         if fix_source_flux[dataset][0] == 0. and fix_source_flux[dataset][1] == 0.:
        #             fix_blend_flux[dataset] = 0.
        return fix_blend_flux

    def _get_flux_from_mag_safe(self, mag):
        """return flux from magnitude, handling NaN values safely."""
        if np.isnan(mag):
            print("Warning: magnitude is NaN, setting flux to 0.")
            return np.float64(0.)
        return mm.Utils.get_flux_from_mag(mag)

    def _get_datasets(self):
        """
        construct a list of MulensModel.MulensData objects
        """
        self._ET_fit = []
        super()._get_datasets()

    def _get_1_dataset(self, file_, kwargs):
        """
        Construct a single dataset and possibly rescale uncertainties in it.
        """
        scaling = file_.pop("scale_errorbars", None)
        _ = file_.pop("fit_errorbars", None)
        bad = file_.pop("bad", None)
        ET_fit = file_.pop("ET_fit", False)
        _ = file_.pop("MIST_bandpass", None)
        try:
            dataset = mm.MulensData(**{**kwargs, **file_})
        except FileNotFoundError:
            raise FileNotFoundError(
                'Provided file path does not exist: ' +
                str(file_['file_name']))
        except Exception:
            print('Something went wrong while reading file ' +
                  str(file_['file_name']), file=sys.stderr)
            raise
        if scaling is not None:
            dataset.scale_errorbars(**scaling)
        if bad is not None:
            self._parse_bad(bad, dataset)

        if dataset.ephemerides_file is not None:
            self._satellites_names.append(dataset.plot_properties['label'])
            self._satellites_colors.append(dataset.plot_properties['color'])

        if ET_fit is False:
            self._ET_fit.append(False)
        else:
            self._ET_fit.append(True)

        return dataset

    def _print_yaml_best_model(self, begin="", mode=None):
        super()._print_yaml_best_model()
        self._print_yaml_track(begin=begin, mode=mode)

    def _print_yaml_track(self, begin="", mode=None):
        yaml_txt = begin + "Sources parameters:\n"

        parameters_star = self._other_parameters_dict
        if self._fixed_parameters is not None:
            parameters_star = {**parameters_star, **self._fixed_parameters}

        (mags, self._source_parameters, failed) = self._get_seds(self._best_model_theta, parameters_star)

        for (i, parameters) in enumerate(self._source_parameters):
            yaml_txt += (begin + "  Source_{:d}:\n").format(i+1)
            for key, value in parameters.items():
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                yaml_txt += (begin + "    {:s}: {:s}\n").format(key, str(value))
            for band, mag in mags.items():
                yaml_txt += (begin + "    {:s}_{:s}: {:s}\n").format(band, 'mag', str(mag[i]))

        yaml_txt += begin + "chi2 per dataset:\n"
        for (i, dataset) in enumerate(self._datasets):
            chi2 = self._event.get_chi2_for_dataset(i)
            yaml_txt += (begin + "  {:s}: {:f}\n").format(dataset.plot_properties['label'], chi2)
        print(yaml_txt, end="", **self._yaml_kwargs)

    def _return_ln_prob(self, value, fluxes=None):
        """
        used to parse output of _ln_prob() in order to make that function
        shorter
        """
        if value == -np.inf:
            n_source_params = self._model.n_sources * (len(self._parameters_star_aux) + 2)
            if self._reparametrized_MM_parameters_values is not None:
                if self._return_fluxes:
                    return (value, [0.] * self._n_fluxes + [0.] * len(self._reparametrized_MM_parameters) +
                            [0.] * n_source_params)
                else:
                    return (value, [0.] * len(self._reparametrized_MM_parameters) + [0.] * n_source_params)
            else:
                if self._return_fluxes:
                    return (value, [0.] * self._n_fluxes + [0.] * n_source_params)
                else:
                    return (value, [0.] * n_source_params)
        else:
            sources_values = list(self.sources_values)
            if self._return_fluxes:
                if fluxes is None:
                    raise ValueError('Unexpected error!')
                if self._reparametrized_MM_parameters_values is not None:
                    return (value, fluxes + self._reparametrized_MM_parameters_values + sources_values)
                else:
                    return (value, fluxes + sources_values)
            else:
                if self._reparametrized_MM_parameters_values is not None:
                    return (value, self._reparametrized_MM_parameters_values + sources_values)
                else:
                    return (value, sources_values)

    def _get_fluxes_to_print_EMCEE(self):
        """
        prepare values to be printed for EMCEE fitting
        """
        if self._reparametrized_MM_parameters is None:
            n_reparametrized_MM_parameters = 0
        else:
            n_reparametrized_MM_parameters = len(self._reparametrized_MM_parameters)
        try:
            blob_samples = np.array(self._sampler.get_blobs(flat=True, discard=self._fitting_parameters['n_burn']))
        except Exception as exception:
            raise ValueError('There was some issue with blobs:\n' +
                             str(exception))
        self._blob_reparametrized = blob_samples[:, self._n_fluxes:self._n_fluxes + n_reparametrized_MM_parameters]
        self._blob_sources_values = blob_samples[:, self._n_fluxes + n_reparametrized_MM_parameters:]
        blob_samples = blob_samples[:, :self._n_fluxes]

        return blob_samples

    def _generate_random_parameters(self):
        """
        Generate a number of starting parameters values.
        It is checked if parameters are within the prior.
        """
        max_iteration = 20 * self._n_walkers
        if self._fit_constraints["no_negative_blending_flux"]:
            max_iteration *= 5
        max_iteration *= 5  # extra factor for EPP
        starting = []
        for parameter in self._fit_parameters:
            settings = self._starting_parameters[parameter]
            values = self._get_samples_from_distribution(
                max_iteration, settings)
            starting.append(values)

        return np.array(starting).T.tolist()

    def _parse_results(self):
        """
        Call the function that prints and saves results
        """
        if self._fit_method == "EMCEE":
            self._parse_results_EMCEE()
            if self._posterior_file_name is not None:
                self._save_posterior_EMCEE()
        elif self._fit_method == "MultiNest":
            self._parse_results_MultiNest()
        elif self._fit_method == "UltraNest":
            self._parse_results_UltraNest()
        else:
            raise ValueError('internal bug')
        if self._reparametrized_MM_parameters is not None:
            self._parse_reparametrized()
        if self._blob_sources_values is not None:
            self._parse_sources()
        # Below we close open files and remove temporary ones.
        if self._yaml_results:
            if self._yaml_results_file is not sys.stdout:
                self._yaml_results_file.close()
        if self._fit_method == "MultiNest":
            if self._MN_temporary_files:
                shutil.rmtree(self._kwargs_MultiNest['outputfiles_basename'],
                              ignore_errors=True)

    def _parse_sources(self, mode=None):
        """
        Printing results for source parameters from isochrone models
        """
        if self.sources_dict is not None:
            ids = []
            for key_source in self.sources_dict.keys():
                for key_param in self.sources_dict[key_source].keys():
                    ids.append(f'{key_source}_{key_param}')
            data = self._blob_sources_values

            if self._fit_method == "EMCEE":
                results = self._get_weighted_percentile(data)
            elif self._fit_method in ["MultiNest", "UltraNest"]:
                if mode is None:
                    weights = self._samples_flat_weights
                else:
                    weights = self._samples_modes_flat_weights[mode]
                results = self._get_weighted_percentile(data, weights=weights)

            begin = "    "
            print("Sources:", **self._yaml_kwargs)
            print("# [median, sigma+, sigma-]", **self._yaml_kwargs)
            print(self._format_results(ids, results, yaml=True, begin=begin), **self._yaml_kwargs)

            begin = "    "
            print("Sources:")
            print("# [median, sigma+, sigma-]")
            print(self._format_results(ids, results))

    def _plot_track(self):
        """
        Plot the track of the source in the Hertzsprung-Russell diagram.
        """
        if self._plot_track_file is not None:
            if len(self._MIST_bandpass) < 2:
                print('Not enough bandpasses to plot the track.')
                print('Setting to I and V bandpasses for plotting.')
                self._setup_brutus(bandpass=['Bessell_I', 'Bessell_V'])

            plt.figure(figsize=(10, 10))
            plt.xlabel('{:s} - {:s}'.format(self._MIST_bandpass[1], self._MIST_bandpass[0]))
            plt.ylabel('{:s}'.format(self._MIST_bandpass[0]))

            parameters_star = self._other_parameters_dict
            if self._fixed_parameters is not None:
                parameters_star = {**parameters_star, **self._fixed_parameters}

            (mags, self._source_parameters, failed) = self._get_seds(self._best_model_theta, parameters_star)

            colors = ['red', 'orange']
            for (i, parameters) in enumerate(self._source_parameters):
                plt.scatter(mags[self._MIST_bandpass[1]][i] - mags[self._MIST_bandpass[0]][i],
                            mags[self._MIST_bandpass[0]][i], label=f'Source {i+1}', zorder=3,
                            color=colors[i], alpha=0.8)
            self.plot_whole_track(parameters_star)
            self.plot_neighbourhood()
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.colorbar(label='EEP')
            plt.title('EEP Track')
            plt.legend()
            file = self._plot_track_file
            plt.savefig(file, dpi=300)

    def plot_whole_track(self, parameters_star):
        """
        Plot the whole evolutionary track on the CMD
        """
        eep_grid = np.linspace(202, 808, 500)
        _all = np.full(
            (self._model.n_sources, 3, len(eep_grid)), np.nan)
        for i, eep in enumerate(eep_grid):
            (mags, source_parameters, failed) = self._get_seds(self._best_model_theta,
                                                               parameters_star, EEP_S=eep, loga_max=20.)
            for j in range(self._model.n_sources):
                if not np.isnan(mags[self._MIST_bandpass[1]][j]) and not np.isnan(mags[self._MIST_bandpass[0]][j]):
                    _all[j][0][i] = mags[self._MIST_bandpass[1]][j] - mags[self._MIST_bandpass[0]][j]
                    _all[j][1][i] = mags[self._MIST_bandpass[0]][j]
                    _all[j][2][i] = eep
        for j in range(self._model.n_sources):
            plt.scatter(_all[j][0], _all[j][1], c=_all[j][2], zorder=2, alpha=0.8)

    def plot_neighbourhood(self):
        """
        Plot the field stars in the CMD if a photometric map file is provided.
        """
        if self._photometric_map_file is None:
            print("""<photometric_map_file> is not set
The color-magnitude diagram (CMD) will include only the evolutionary tracks and will not show field stars.
To include field stars in the CMD, set photometric_map_file to a file containing field-star magnitudes
in the same photometric bands used by the MIST models.""")
            return

        data = np.genfromtxt(self._photometric_map_file, delimiter=' ', names=True)
        plt.scatter(data[self._MIST_bandpass[1]] - data[self._MIST_bandpass[0]],
                    data[self._MIST_bandpass[0]], alpha=0.5, color='gray', s=0.1)

    def run_fit(self):
        """
        Run the fit, print the output, and make the plots.

        This function does not accept any parameters. All the settings
        are passed via __init__().
        """
        if self._task != "fit":
            raise ValueError('wrong settings to run .run_fit()')

        self._check_plots_parameters()
        self._check_model_parameters()
        self._check_other_fit_parameters()
        self._parse_other_output_parameters()
        self._setup_brutus()
        self._get_datasets()
        self._check_ulens_model_parameters()
        self._get_parameters_ordered()
        self._parse_errorbars_fit_params()
        self._get_parameters_latex()
        self._set_prior_limits()

        if self._fit_method == "EMCEE":
            self._parse_starting_parameters()

        self._check_fixed_parameters()
        # self._parse_degeneracy()  #if AD used
        self._make_model_and_event()
        self._parse_fitting_parameters()
        self._parse_fit_constraints()

        if self._fit_method == "EMCEE":
            self._get_starting_parameters()

        self._setup_fit()
        self._run_fit()
        self._finish_fit()
        self._parse_results()
        self._write_residuals()
        self._make_plots()
        self._plot_track()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')
    if 'yaml' in import_failed:
        raise ImportError('module "yaml" could not be imported :(')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = UlensModelFitEvolutionaryTracks(**settings)
    self = ulens_model_fit
    ulens_model_fit.run_fit()
