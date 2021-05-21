import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ColorConverter

from astropy import units as u
from astropy.coordinates import SkyCoord

from MulensModel.modelparameters import ModelParameters
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.trajectory import Trajectory
from MulensModel.caustics import Caustics
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs
from MulensModel.coordinates import Coordinates


class Model(object):
    """
    A Model for a microlensing event with the specified parameters.

    Arguments :
        parameters: *dictionary*,
        :py:class:`~MulensModel.modelparameters.ModelParameters`

            see
            :py:class:`MulensModel.modelparameters.ModelParameters`

        :py:obj:`coords`: *str*, *astropy.SkyCoords*,
        *MulensModel.Coordsinates*, optional

            Sky coordinates of the event. If type is *str*, then it is
            assumed that the units are hour angle and degrees for RA and Dec,
            respectively.

        ra, dec: *str*, optional
            Sky coordinates of the event.

        ephemerides_file: *str*, optional
            Specify name of the file with satellite ephemerides. See
            :py:class:`~MulensModel.mulensdata.MulensData` for more details.
            Note that if you provide file name here, then it will affect
            all calculations for this model. In most cases, you want to
            combine ground-based and satellite data and in those cases set
            ephemerides_file for specific
            :py:class:`~MulensModel.mulensdata.MulensData` instance
            to pass satellite information.

    Attributes :
        ephemerides_file: *str*
            Name of file with satellite ephemerides.

        data_ref: *int* or :py:class:`~MulensModel.mulensdata.MulensData`
            Reference dataset. If *int* then gives index of reference dataset
            in :py:attr:`~datasets`.

    Default values for parallax are all True. Use :py:func:`parallax()`
    to turn different parallax effects ON/OFF. If using satellite
    parallax, you may also specify an `ephemerides_file` (see
    :py:class:`~MulensModel.mulensdata.MulensData`).

    Note that you can print an instance of Model, which shows you parameters
    in a nice way, e.g.,

    .. code-block:: python

      model = Model(parameters={'t_0': 2456789.0, ....})
      print(model)

    """

    def __init__(
            self, parameters=None, coords=None, ra=None, dec=None,
            ephemerides_file=None):

        # Initialize the parameters of the model
        if isinstance(parameters, ModelParameters):
            self._parameters = parameters
        else:
            self._parameters = ModelParameters(parameters)

        # Set the coordinates of the event
        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
        if coords is not None:
            self._coords = Coordinates(coords)

        if ra is not None:
            if dec is not None:
                self._coords = Coordinates(ra, dec)
            else:
                raise AttributeError(coords_msg)
        else:
            if dec is not None:
                raise AttributeError(coords_msg)

        self.ephemerides_file = ephemerides_file
        self._satellite_skycoord = None

        # Set some defaults
        self._parallax = {'earth_orbital': True,
                          'satellite': True,
                          'topocentric': True}
        self._default_magnification_method = 'point_source'
        self._methods = None
        self._methods_parameters = {}
        self._caustics = None
        self.data_ref = None
        self._fit = None

        self._limb_darkening_coeffs = LimbDarkeningCoeffs()
        self._bandpasses = []

        self._source_flux_ratio_constraint = dict()

        self._datasets = []
        self._datasets_set_single = False

    def __repr__(self):
        return '{0}'.format(self.parameters)

    def plot_magnification(
            self, times=None, t_range=None, t_start=None, t_stop=None, dt=None,
            n_epochs=None, subtract_2450000=False, subtract_2460000=False,
            satellite_skycoord=None, gamma=0., flux_ratio_constraint=None,
            **kwargs):
        """
        Plot the model magnification curve.

        Keywords :
            see :py:func:`plot_lc()`

            gamma:
                see :py:func:`magnification()`

            satellite_skycoord:
                see :py:func:`plot_trajectory()`

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.

        """
        if 'fit_blending' in kwargs:
            raise ValueError(
                'Keyword fit_blending is allowed in Model.plot_lc() ' +
                'but not Model.plot_magnification()')
        if times is None:
            times = self.set_times(
                t_range=t_range, t_start=t_start, t_stop=t_stop, dt=dt,
                n_epochs=n_epochs)
        subtract = self._subtract(subtract_2450000, subtract_2460000)

        if satellite_skycoord is not None:
            if isinstance(satellite_skycoord, SatelliteSkyCoord):
                satellite = satellite_skycoord.get_satellite_coords(times)
            elif not isinstance(satellite_skycoord, SkyCoord):
                raise TypeError('Wrong type of satellite_skycoord in ' +
                                'Model.plot_magnification()')
        else:
            satellite = None
        if self.n_sources == 2:
            if (flux_ratio_constraint is None and
                    None not in self._source_flux_ratio_constraint):
                flux_ratio_constraint = (
                    self._flux_ratio_constraint_for_plotting())

        magnification = self.magnification(
            times, satellite_skycoord=satellite, gamma=gamma,
            flux_ratio_constraint=flux_ratio_constraint)

        self._plt_plot(times-subtract, magnification, kwargs)
        plt.ylabel('Magnification')
        plt.xlabel(self._subtract_xlabel(subtract_2450000, subtract_2460000))

    def plot_lc(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, data_ref=None, f_source=None, f_blend=None,
            subtract_2450000=False, subtract_2460000=False,
            flux_ratio_constraint=None, fit_blending=None,
            gamma=None, bandpass=None, **kwargs):
        """
        Plot the model light curve in magnitudes.

        Keywords :
            times: [*float*, *list*, *numpy.ndarray*]
                a list of times at which to plot the magnifications

            t_range, t_start, t_stop, dt, n_epochs: see :py:func:`set_times`

            data_ref: *int* or a
            :py:class:`~MulensModel.mulensdata.MulensData` object

                Reference dataset to scale the model to. See
                :py:func:`get_ref_fluxes()`

            f_source, f_blend: *float*
                Explicitly specify the source and blend fluxes in a
                system where flux = 1 corresponds to
                :obj:`MulensModel.utils.MAG_ZEROPOINT` (= 22 mag).

            subtract_2450000, subtract_2460000: *boolean*, optional
                If True, subtracts 2450000 or 2460000 from the time
                axis to get more human-scale numbers. If using, make
                sure to also set the same settings for all other
                plotting calls (e.g. :py:func:`plot_data()`)

            flux_ratio_constraint: instance of
            :py:class:`~MulensModel.mulensdata.MulensData` or *str*, optional
                Option for binary source models only.
                Data or bandpass to constrain the flux ratio of sources.

            fit_blending: *boolean*
                *True* if blending flux is going to be fitted (default),
                *False* if blending flux is fixed at 0.

            gamma:
                see :py:func:`magnification()`

            bandpass: *str*
                bandpass for defining the limb-darkening coefficient gamma.
                Requires that the limb_darkenging coefficients have been set
                using :py:func:`set_limb_coeff_u()` or
                :py:func:`set_limb_coeff_gamma()`.

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.

        Provide `data_ref` or (`f_source`, `f_blend`) if you want to
        plot in flux units different than last value of `data_ref`
        (defaults to the first dataset).

        """
        if data_ref is not None:
            self.data_ref = data_ref
        if f_source is None and f_blend is None:
            if self.data_ref is None:
                raise ValueError('No reference dataset of fluxes provided. ' +
                                 "If you don't have a dataset, then try " +
                                 "plot_magnification() instead of plot_lc().")
        elif f_source is None or f_blend is None:
            raise AttributeError(
                'If f_source is set, f_blend must also be set and vice versa.')
        else:
            if fit_blending is not None:
                raise ValueError(
                    "Error in plot_lc(). Providing f_source, " +
                    "f_blend, and fit_blending doesn't make sense.")

        if f_source is None and f_blend is None:
            (f_source, f_blend) = self.get_ref_fluxes(
                data_ref=self.data_ref, fit_blending=fit_blending)

        if times is None:
            times = self.set_times(
                t_range=t_range, t_start=t_start, t_stop=t_stop, dt=dt,
                n_epochs=n_epochs)
        elif isinstance(times, list):
            times = np.array(times)

        if (bandpass is not None) and (gamma is not None):
            raise ValueError('Only one of bandpass and gamma can be set')
        elif (bandpass is None) and (gamma is None):
            gamma = 0.
        elif bandpass is not None:
            if bandpass not in self._bandpasses:
                raise KeyError(
                    'No limb-darkening coefficient set for {0}'.format(
                        bandpass))
            else:
                gamma = self.get_limb_coeff_gamma(bandpass)

        else:
            pass

        if self.n_sources == 2:
            if (flux_ratio_constraint is None and
                    None not in self._source_flux_ratio_constraint):
                flux_ratio_constraint = (
                    self._flux_ratio_constraint_for_plotting())

        magnification = self.magnification(
            times, gamma=gamma,
            flux_ratio_constraint=flux_ratio_constraint)
        flux = f_source * magnification + f_blend
        subtract = self._subtract(subtract_2450000, subtract_2460000)

        self._plt_plot(times-subtract, Utils.get_mag_from_flux(flux), kwargs)
        plt.ylabel('Magnitude')
        plt.xlabel(self._subtract_xlabel(subtract_2450000, subtract_2460000))

        (ymin, ymax) = plt.gca().get_ylim()
        if ymax > ymin:
            plt.gca().invert_yaxis()

    def _flux_ratio_constraint_for_plotting(self):
        """
        Warning or error message for lack of source flux ratio information.
        """
        msg = (
                '\n\nTo plot magnification for binary source model you have ' +
                'to set the flux ratio (using set_source_flux_ratio()) ' +
                'or provide dataset which will be used to find flux ' +
                'ratio (keyword flux_ratio_constraint).\n\n')
        if len(self._datasets) == 1:
            warnings.warn(msg + 'You have provided only one dataset, so ' +
                          "for now we will use that one.", RuntimeWarning)
            return self._datasets[0]
        elif len(self._datasets) > 1:
            warnings.warn(msg + 'You have provided only some datasets and ' +
                          "for now we will use the first one.", RuntimeWarning)
            return self._datasets[0]
        else:
            raise ValueError(
                'Not enough information to plot the model magnification. ' +
                'Use set_source_flux_ratio() function ' +
                'or flux_ratio_constraint keyword.')

    def _plt_plot(self, x, y, kwargs):
        """
        safe run of matplotlib.pyplot.plot()
        """
        try:
            plt.plot(x, y, **kwargs)
        except Exception:
            print("kwargs passed to plt.plot():")
            print(kwargs)
            raise

    def _subtract(self, subtract_2450000, subtract_2460000):
        """
        find value of HJD to be subtracted
        """
        if subtract_2450000:
            if subtract_2460000:
                raise ValueError('subtract_2450000 and subtract_2460000 ' +
                                 'cannot be both True')
            subtract = 2450000.
        elif subtract_2460000:
            subtract = 2460000.
        else:
            subtract = 0.
        return subtract

    def _subtract_xlabel(self, subtract_2450000, subtract_2460000):
        """
        string that would be past to plt.xlabel()
        """
        if subtract_2450000:
            if subtract_2460000:
                raise ValueError('subtract_2450000 and subtract_2460000 ' +
                                 'cannot be both True')
            out = 'Time - 2450000'
        elif subtract_2460000:
            out = 'Time - 2460000'
        else:
            out = 'Time'
        return out

    def reset_plot_properties(self):
        """
        This function will be **deprecated**.

        Resets internal plotting properties of all attached datasets.
        """
        warnings.warn('reset_plot_properties() will be deprecated in future',
                      FutureWarning)
        for data in self._datasets:
            data.plot_properties = {}

    def _set_default_colors(self):
        """
        If the user has not specified a color for a dataset, assign
        one.
        """
        colors = [cycle['color'] for cycle in rcParams['axes.prop_cycle']]

        # Below we change the order of colors to most distinct first.
        used_colors = []
        for data in self._datasets:
            if 'color' in data.plot_properties.keys():
                used_colors.append(data.plot_properties['color'])
        if len(used_colors) == len(self._datasets):
            return
        if len(used_colors) == 0:
            differences = None
        else:
            d_col = self._color_differences
            diffs = np.array([np.min(d_col(used_colors, c)) for c in colors])
            indexes = np.argsort(diffs)[::-1]
            colors = [colors[i] for i in indexes]
            differences = diffs[indexes]

        # Assign colors when needed.
        color_index = 0
        for data in self._datasets:
            if 'color' not in data.plot_properties.keys():
                if differences is not None:
                    if differences[color_index] < 0.35:
                        msg = ('The color assign to one of the datasets in ' +
                               'automated way (' + colors[color_index] +
                               ') is very similar to already used color')
                        warnings.warn(msg, UserWarning)
                data.plot_properties['color'] = colors[color_index]
                color_index += 1
                if color_index == len(colors):
                    color_index = 0
                    msg = ('Too many datasets without colors assigned - ' +
                           'same color will be used for different datasets')
                    warnings.warn(msg, UserWarning)

    def _color_differences(self, color_list, color):
        """
        Calculate color difference between a list of colors and a single color.
        Uses algorithm from
        `this Wikipedia page<https://en.wikipedia.org/wiki/Color_difference>`_.
        Arguments :
            color_list: *list* of *str*
                list of matplotlib colors e.g., ``['black', '#E9967A']``
            color: *str*
                single matplotlib color
        Returns :
            differences: *np.ndarray*
                differences of colors, values < 0.3 are very similar
        """
        rgba = ColorConverter.to_rgba
        array = np.array(
            [[float(x) for x in list(rgba(c))[:3]] for c in color_list])
        # We use float above because some versions of matplotlib return str.
        color_value = [float(x) for x in list(rgba(color))[:3]]
        mean_red = 0.5 * (array[:, 0] + color_value[0])
        diffs = (array - color_value)**2
        add_1 = (2. + mean_red) * diffs[:, 0]
        add_2 = 4. * diffs[:, 1]
        add_3 = (3. + mean_red) * diffs[:, 2]
        return np.sqrt(add_1 + add_2 + add_3)

    def _check_old_plot_kwargs(self, **kwargs):
        """
        Check for deprecated "_list" keywords. Issue a warning, then
        transfer the properties to the new
        :py:attr:`mulensdata.MulensData.plot_properties` system.
        """
        old_plot_keywords = [
            'color_list', 'marker_list', 'size_list',
            'label_list', 'alpha_list', 'zorder_list']

        for old_keyword in old_plot_keywords:
            if kwargs[old_keyword] is not None:
                warnings.warn('Keyword "' + old_keyword + '" is deprecated.' +
                              ' Use MulensData.plot_properties instead.',
                              FutureWarning)
                values = kwargs[old_keyword]
                key = old_keyword[:-5]
                for (dataset, value) in zip(self._datasets, values):
                    dataset.plot_properties[key] = value

    def plot_data(
            self, data_ref=None, show_errorbars=None, show_bad=None,
            color_list=None, marker_list=None, size_list=None,
            label_list=None, alpha_list=None, zorder_list=None,
            subtract_2450000=False, subtract_2460000=False, **kwargs):
        """
        Plot the data scaled to the model.

        Keywords (all optional):
            data_ref: see :py:func:`get_ref_fluxes()`
                If data_ref is not specified, uses the first dataset
                as the reference for flux scale.

            show_errorbars: *boolean* or *None*
                Do you want errorbars to be shown for all datasets?
                Default is *None*, which means the option is taken from each
                dataset plotting properties (for which default is *True*).
                If *True*, then data are plotted using matplotlib.errorbar().
                If *False*, then data are plotted using matplotlib.scatter().

            show_bad: *boolean* or *None*
                Do you want data marked as bad to be shown?
                Default is *None*, which means the option is taken from each
                dataset plotting properties (for which default is *False*).
                If bad data are shown, then they are plotted with 'x' marker.

            subtract_2450000, subtract_2460000: *boolean*
                If True, subtracts 2450000 or 2460000 from the time
                axis to get more human-scale numbers. If using, make
                sure to also set the same settings for all other
                plotting calls (e.g. :py:func:`plot_lc()`).

            ``**kwargs``:
                Passed to matplotlib plotting functions. Contrary to
                previous behavior, ``**kwargs`` are no longer remembered.

        """

        self._check_old_plot_kwargs(
            color_list=color_list, marker_list=marker_list,
            size_list=size_list, label_list=label_list,
            alpha_list=alpha_list, zorder_list=zorder_list)
        self._set_default_colors()

        if data_ref is not None:
            self.data_ref = data_ref

        # Set plot limits
        t_min = 3000000.
        t_max = 0.
        subtract = self._subtract(subtract_2450000, subtract_2460000)

        # Get fluxes for all datasets
        fit = Fit(data=self._datasets, magnification=self.data_magnification)
        fit.fit_fluxes()
        self._fit = fit

        for data in self._datasets[::-1]:
            data.plot(
                phot_fmt='mag', show_errorbars=show_errorbars,
                show_bad=show_bad, subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000, model=self,
                **kwargs)

            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))

        # Plot properties
        plt.ylabel('Magnitude')
        plt.xlabel(self._subtract_xlabel(subtract_2450000, subtract_2460000))
        plt.xlim(t_min-subtract, t_max-subtract)

        (ymin, ymax) = plt.gca().get_ylim()
        if ymax > ymin:
            plt.gca().invert_yaxis()

    def plot_residuals(
            self, show_errorbars=None,
            color_list=None, marker_list=None, size_list=None,
            label_list=None, alpha_list=None, zorder_list=None,
            data_ref=None, subtract_2450000=False, subtract_2460000=False,
            show_bad=None, **kwargs):
        """
        Plot the residuals (in magnitudes) of the model.

        For explanation of keywords, see doctrings in
        :py:func:`plot_data()`. Note different order of keywords.
        """

        self._check_old_plot_kwargs(
            color_list=color_list, marker_list=marker_list,
            size_list=size_list, label_list=label_list,
            alpha_list=alpha_list, zorder_list=zorder_list)
        self._set_default_colors()

        if data_ref is not None:
            self.data_ref = data_ref

        # Plot limit parameters
        t_min = 3000000.
        t_max = 0.
        subtract = self._subtract(subtract_2450000, subtract_2460000)

        # Plot zeropoint line
        plt.plot([0., 3000000.], [0., 0.], color='black')

        # Plot residuals
        for data in self._datasets[::-1]:
            data.plot(
                phot_fmt='mag', show_errorbars=show_errorbars,
                show_bad=show_bad, subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000, model=self,
                plot_residuals=True, **kwargs)
            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))

        # Plot properties
        y_lim = np.max([np.abs(y_lim) for y_lim in plt.gca().get_ylim()])
        if y_lim > 1.:
            y_lim = 0.5
        plt.ylim(y_lim, -y_lim)
        plt.xlim(t_min-subtract, t_max-subtract)
        plt.ylabel('Residuals')
        plt.xlabel(self._subtract_xlabel(subtract_2450000, subtract_2460000))

    def get_residuals(self, data_ref=None, type='mag', data=None):
        """
        Calculate the residuals from the model for
        each dataset at once, or just a single dataset.

        Note: if residuals are returned in magnitudes, they are
        transformed to the magnitude system specified by `data_ref`,
        so only suitable for plotting.

        Keywords :
            data_ref: optional
                see :py:func:`get_ref_fluxes()`

            type: *str*, optional
                specify whether the residuals should be returned in
                magnitudes ('mag') or in flux ('flux'). Default is
                'mag'.

            data: :py:class:`~MulensModel.mulensdata.MulensData`, optional
                dataset for which residuals are returned. If specified,
                then returned lists are single element.

        Returns :
            residuals: *list*
                each element of the list is a *np.ndarray* with the
                residuals for the corresponding dataset.

            errorbars: *list*
                the scaled errorbars for each point. For plotting
                errorbars for the residuals.
        """
        if data_ref is not None:
            self.data_ref = data_ref
        data_ref_ = self._get_data_ref(data_ref)

        if data is not None:
            data_list = [data]
            fit_data = list(set([data, data_ref_]))
            magnifications = [self.get_data_magnification(d) for d in fit_data]
            fit = Fit(data=fit_data, magnification=magnifications)
        else:
            data_list = self._datasets
            fit = Fit(data=data_list, magnification=self.data_magnification)
        fit.fit_fluxes()
        self._fit = fit

        residuals = []
        errorbars = []
        # Calculate residuals.
        for data_ in data_list:
            f_source = self.fit.flux_of_sources(data_)
            f_blend = self.fit.blending_flux(data_)
            if type == 'mag':
                f_source_0 = self.fit.flux_of_sources(data_ref_)
                f_blend_0 = self.fit.blending_flux(data_ref_)
                magnification = self.get_data_magnification(data_)
                model_mag = Utils.get_mag_from_flux(
                    f_blend_0 + f_source_0 * magnification)
                flux = (f_source_0 * (data_.flux - f_blend) /
                        f_source + f_blend_0)
                err_flux = f_source_0 * data_.err_flux / f_source
                (mag, err) = Utils.get_mag_and_err_from_flux(flux, err_flux)
                residuals.append(mag - model_mag)
                errorbars.append(err)
            elif type == 'flux':
                magnification = self.get_data_magnification(data_)
                model_flux = f_blend + f_source * magnification
                residuals.append(data_.flux - model_flux)
                errorbars.append(data_.err_flux)
            else:
                raise ValueError("type keyword must be either 'mag' or 'flux'")

        return (residuals, errorbars)

    def plot_caustics(self, n_points=5000, epoch=None, **kwargs):
        """
        Plot the caustic structure. See
        :py:func:`MulensModel.caustics.Caustics.plot()` for binary lenses.
        For a single lens it just marks (0, 0) point and
        the first two parameters are ignored.

        Additional parameters :
            n_points: *int*, optional
                The number of points to calculate along the caustic.
                Defaults to 5000.

            epoch: *float*, optional
                Epoch for which separation *s* will be used. Important
                for models with orbital motion. Defaults to *t_0_kep*,
                which defaults to *t_0*.

            ``**kwargs``:
                keywords accepted by :py:func:`matplotlib.pyplot.scatter()`
        """
        if self.n_lenses == 1:
            plt.scatter([0], [0], **kwargs)
        else:
            self.update_caustics(epoch=epoch)
            self.caustics.plot(n_points=n_points, **kwargs)

    def update_caustics(self, epoch=None):
        """
        Updates :py:attr:`~caustics` property for given epoch.

        Parameters :
            epoch: *float*
                For orbital motion models, epoch for which separation *s*
                is calculated to calculate :py:attr:`~caustics`. Defaults
                to *t_0_kep*, which defaults to *t_0*.
        """
        if epoch is None:
            s = self.parameters.s
        else:
            s = self.parameters.get_s(epoch)

        if self._caustics is not None:
            if s == self._caustics.s and self.parameters.q == self._caustics.q:
                return

        self._caustics = Caustics(q=self.parameters.q, s=s)

    def plot_trajectory(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, caustics=False, show_data=False,
            arrow=True, satellite_skycoord=None, arrow_kwargs=None, **kwargs):
        """
        Plot the source trajectory.

        Keywords (all optional) :

            times, t_range, t_start, t_stop, dt, n_epochs:
                May all be used to specify exactly when to plot the
                source trajectory. See also :py:func:`plot_lc()` and
                :py:func:`set_times()`.

            caustics: *boolean*
                plot the caustic structure in addition to the source
                trajectory. default=False (off). For finer control of
                plotting features, e.g. color, use :py:func:`plot_caustics()`
                instead.

            show_data: *boolean*
                Mark epochs of data. Use :py:func:`plot_source_for_datasets()`
                for more control of how data epochs are marked.

            arrow: *boolean*
                Show the direction of the source motion. Default is *True*.

            satellite_skycoord: *astropy.SkyCoord* or
            :py:class:`MulensModel.satelliteskycoord.SatelliteSkyCoord`

                Allows the user to specify that the trajectory is calculated
                for a satellite. If *astropy.SkyCoord* object is provided,
                then these are satellite positions for all epochs.
                See also :py:func:`get_satellite_coords()`

            arrow_kwargs: *dict*
                Kwargs that are passed to :py:func:`pyplot.arrow()`. If no
                color is given here, then we use one specified in ``**kwargs``
                and if nothing is there, then we use black. The size of
                the arrow is determined based on limits of current axis.
                If those are not adequate, then change the size by specifying
                *width* keyword and maybe other as well. Note that
                *arrow_kwargs* are of *dict* type and are different than
                ``**kwargs``.

            ``**kwargs``
                Controls plotting features of the trajectory. It's passed to
                :py:func:`pyplot.plot()`.

        Note that in order to have equal scaling of both axis
        (i.e., make circles look circular), you have to call appropriate
        *pyplot* command. This can be one of these commands:

        .. code-block:: python

          pyplot.axis('equal')
          pyplot.axis('scaled')
          pyplot.axis('square')
          pyplot.gca().set_aspect('equal')

        They have slightly different behavior.

        """
        if not arrow and arrow_kwargs is not None:
            raise ValueError(
                "arrow_kwargs can be only given if arrow is True")

        if times is None:
            times = self.set_times(
                t_range=t_range, t_start=t_start, t_stop=t_stop, dt=dt,
                n_epochs=n_epochs)
        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(times)
        else:
            if isinstance(satellite_skycoord, SatelliteSkyCoord):
                satellite_skycoord = satellite_skycoord.get_satellite_coords(
                    times)
            elif not isinstance(satellite_skycoord, SkyCoord):
                raise TypeError('Wrong type of satellite_skycoord in ' +
                                'Model.plot_trajectory()')

        if self.n_sources == 1:
            self._plot_single_trajectory(
                times, self.parameters, satellite_skycoord,
                arrow, arrow_kwargs, **kwargs)
        elif self.n_sources == 2:
            self._plot_single_trajectory(
                times, self.parameters.source_1_parameters,
                satellite_skycoord, arrow, arrow_kwargs, **kwargs)
            self._plot_single_trajectory(
                times, self.parameters.source_2_parameters,
                satellite_skycoord, arrow, arrow_kwargs, **kwargs)
        else:
            raise ValueError(
                    'Wrong number of sources: {:}'.format(self.n_sources))

        if caustics:
            self.plot_caustics(marker='.', color='red')

        if show_data:
            self.plot_source_for_datasets()

    def _plot_single_trajectory(self, times, parameters, satellite_skycoord,
                                arrow, arrow_kwargs, **kwargs):
        """
        Plots trajectory of a single source.
        """
        if len(times) < 2:
            raise ValueError('Trajectory can be plotted when at least two ' +
                             'epochs are provided.')
        trajectory = Trajectory(
            times, parameters=parameters, parallax=self._parallax,
            coords=self._coords, satellite_skycoord=satellite_skycoord)

        self._plt_plot(trajectory.x, trajectory.y, kwargs)

        if arrow:
            if len(times) > 2:
                index = int(len(times)/2)
            else:
                index = 0
            x_0 = trajectory.x[index]
            y_0 = trajectory.y[index]
            d_x = trajectory.x[index+1] - x_0
            d_y = trajectory.y[index+1] - y_0
            dd = 1e6 * (d_x*d_x + d_y*d_y)**.5

            xlim = plt.xlim()
            ylim = plt.ylim()
            width = np.abs(xlim[1]-xlim[0]) * np.abs(ylim[1]-ylim[0])
            width = width**.5 / 100.

            color = kwargs.get('color', 'black')
            kwargs_ = {'width': width, 'color': color, 'lw': 0,
                       'zorder': -np.inf}
            if arrow_kwargs is not None:
                kwargs_.update(arrow_kwargs)
            plt.arrow(x_0, y_0, d_x/dd, d_y/dd, **kwargs_)

    def plot_source(self, times=None, **kwargs):
        """
        Plot source: circles of the radius rho at positions corresponding to
        source positions at times. When the rho is not defined, then X symbols
        are plotted.

        Parameters:
            times: *float* or *np.ndarray*
                epochs for which source positions will be plotted

            ``**kwargs``:
                Keyword arguments passed to matplotlib.Circle_.
                Examples: ``color='red'``, ``fill=False``,
                ``linewidth=3``, ``alpha=0.5``. When the rho is not defined,
                then keyword arguments are passed to matplotlib.plot_.

        Note that it is likely that with default axis scaling, the circles may
        be plotted as ellipses. To mitigate it, use:
        ``plt.gca().set_aspect('equal')`` or ``plt.axis('equal')``
        (the other possible options are ``'scaled'`` or ``'square'``).

        .. _matplotlib.Circle:
          https://matplotlib.org/api/_as_gen/matplotlib.patches.Circle.html
        .. _matplotlib.plot:
          https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html
        """
        if isinstance(times, float):
            times = [times]
        if isinstance(times, list):
            times = np.array(times)

        kwargs_ = {
            'times': times, 'parallax': self._parallax, 'coords': self._coords,
            'satellite_skycoord': self.get_satellite_coords(times)}

        if self.n_sources == 1:
            trajectory = Trajectory(parameters=self.parameters, **kwargs_)
            self._plot_source_for_trajectory(trajectory, **kwargs)
        elif self.n_sources == 2:
            trajectory = Trajectory(
                parameters=self.parameters.source_1_parameters, **kwargs_)
            self._plot_source_for_trajectory(trajectory, **kwargs)
            trajectory = Trajectory(
                parameters=self.parameters.source_2_parameters, **kwargs_)
            self._plot_source_for_trajectory(trajectory, **kwargs)
        else:
            raise ValueError('Wrong number of sources!')

    def plot_source_for_datasets(self, **kwargs):
        """
        Plot source positions for all linked datasets. Colors used for
        each dataset are the same used for plotting photometry.

        Parameters:
            ``**kwargs``:
                see :py:func:`plot_source`
        """
        if 'color' in kwargs:
            raise ValueError(
                'Color cannot be passed to plot_source_for_datasets(). To ' +
                'change color, just change color for given MulensData object')
        if self.datasets is None:
            raise ValueError(
                'No data linked to Model, so calling ' +
                'plot_source_for_datasets() makes no sense.')

        self._set_default_colors()
        if isinstance(self.datasets, MulensData):
            datasets = [self.datasets]
        else:
            datasets = self.datasets

        kwargs_ = dict(kwargs)
        if 'fill' not in kwargs_:
            kwargs_['fill'] = False

        for data in datasets:
            kwargs_['color'] = data.plot_properties['color']
            self.plot_source(times=data.time, **kwargs_)

    def _plot_source_for_trajectory(self, trajectory, **kwargs):
        """
        Internal function for plotting sources.
        """
        axis = plt.gca()

        plot_circle = (trajectory.parameters.rho is not None)
        if not plot_circle:
            warnings.warn(
                'rho is not defined, so X is used to show source positions')

        if 'color' not in kwargs:
            kwargs['color'] = 'red'
        if 'zorder' not in kwargs:
            kwargs['zorder'] = np.inf
        if plot_circle and 'radius' not in kwargs:
            kwargs['radius'] = trajectory.parameters.rho

        xlim = plt.xlim()
        ylim = plt.ylim()
        if plot_circle and not (xlim == (0., 1.) and ylim == (0., 1.)):
            # We've checked that default limits were changed.
            limits = max(np.abs(xlim[1]-xlim[0]), np.abs(ylim[1]-ylim[0]))
            if limits > 1e3 * kwargs['radius']:
                warnings.warn(
                    "Model.plot_source() - the source radius is much smaller" +
                    " than the axis range (at this point). The plotted " +
                    "source may be hard to see. To correct it you can change" +
                    " axis range (plt.xlim() or ylim()) or plot larger " +
                    "circle by passing *radius* via kwargs.",
                    UserWarning)

        if not plot_circle:
            plt.plot(trajectory.x, trajectory.y, 'x', **kwargs)
        else:
            for (x, y) in zip(trajectory.x, trajectory.y):
                axis.add_artist(plt.Circle((x, y), **kwargs))

    def get_ref_fluxes(self, data_ref=None, fit_blending=None):
        """
        Get source and blending fluxes for the model by finding the
        best-fit values compared to data_ref.

        Parameters :
            data_ref: :py:class:`~MulensModel.mulensdata.MulensData` or *int*
                Reference dataset. If *int*, corresponds to the index of
                the dataset in self.datasets. If None, than the first dataset
                will be used.

            fit_blending: *boolean*
                *True* if blending flux is going to be fitted (default),
                *False* if blending flux is fixed at 0.

        Returns :
            f_source: *np.ndarray*
                Sources' flux; normally of size (1). If it is of size (1)
                for a double source model, then it is a sum of fluxes
                of both sources.
            f_blend: *float*
                blending flux

        Determine the reference flux system from the datasets. The
        *data_ref* may either be a dataset or the index of a dataset
        (if :py:func:`Model.set_datasets()` was previously called). If
        *data_ref* is not set, it will use the first dataset. If you
        call this without calling :py:func:`set_datasets()` first,
        there will be an exception and that's on you.
        """
        data = self._get_data_ref(data_ref)

        mags = self.get_data_magnification(data)
        fit = Fit(data=data, magnification=mags)
        if fit_blending is not None:
            fit.fit_fluxes(fit_blending=fit_blending)
        else:
            fit.fit_fluxes()
        self._fit = fit

        f_source = fit.flux_of_sources(data)
        f_blend = fit.blending_flux(data)

        return (f_source, f_blend)

    def _get_data_ref(self, data_ref):
        """
        Guess which reference dataset is talked about. Returns MulensData
        instance.
        """
        if data_ref is None:
            if len(self._datasets) == 0:
                raise ValueError(
                    'You cannot get reference flux for Model if' +
                    ' you have not linked data first.')
            if isinstance(self.data_ref, MulensData):
                data = self.data_ref
            else:
                data = self._datasets[self.data_ref]
        elif isinstance(data_ref, MulensData):
            data = data_ref
            self.data_ref = data_ref
        else:
            data = self._datasets[data_ref]
            self.data_ref = data_ref

        return data

    def set_times(
            self, t_range=None, t_start=None, t_stop=None, dt=None,
            n_epochs=1000):
        """
        Return a list of times. If no keywords are specified, default
        is 1000 epochs from [:math:`t_0 - 1.5 * t_E`, :math:`t_0 + 1.5 * t_E`]
        range.
        For binary source models, respectively, smaller and larger of
        `t_0_1`/`t_0_2` values are used.

        Parameters (all optional):
            t_range: [*list*, *tuple*]
                A range of times of the form [t_start, t_stop]

            t_start, t_stop: *float*
                a start or stop time.

            dt: *float*
                the interval spacing between successive points

            n_epochs: *int*
                the number of epochs (evenly spaced)

        Returns :
            times: *np.ndarray*
                Vector of epochs.
        """
        if t_range is not None:
            if t_start is not None or t_stop is not None:
                raise ValueError(
                    'Model.set_times() - you cannot set t_range and either ' +
                    't_start or t_stop')
            t_start = t_range[0]
            t_stop = t_range[1]

        n_tE = 1.5
        if t_start is None:
            if self.n_sources == 1:
                t_0 = self.parameters.t_0
            else:
                t_0 = min(self.parameters.source_1_parameters.t_0,
                          self.parameters.source_2_parameters.t_0)
            t_start = t_0 - (n_tE * self.parameters.t_E)
        if t_stop is None:
            if self.n_sources == 1:
                t_0 = self.parameters.t_0
            else:
                t_0 = max(self.parameters.source_1_parameters.t_0,
                          self.parameters.source_2_parameters.t_0)
            t_stop = t_0 + (n_tE * self.parameters.t_E)

        if dt is None:
            if n_epochs is None:
                n_epochs = 1000
            n_epochs -= 1
            dt = (t_stop - t_start) / float(n_epochs)

        out = np.arange(t_start, t_stop+dt, dt)
        if out[-1] > t_stop:  # This may happen due to rounding errors.
            out = out[:-1]
        return out

    def set_magnification_methods(self, methods, source=None):
        """
        Sets methods used for magnification calculation. See
        :py:class:`~MulensModel.magnificationcurve.MagnificationCurve`
        for a list of implemented methods.

        Parameters :
            methods: *list*
                List that specifies which methods (*str*) should be
                used when (*float* values for Julian dates). Given
                method will be used for times between the times
                between which it is on the list, e.g.,

                .. code-block:: python

                  methods = [
                      2455746., 'Quadrupole', 2455746.6, 'Hexadecapole',
                      2455746.7, 'VBBL', 2455747., 'Hexadecapole',
                      2455747.15, 'Quadrupole', 2455748.]

            source: *int* or *None*
                Which source given methods apply to? Accepts 1, 2, or *None*
                (i.e., all sources).
        """
        if not isinstance(methods, list):
            raise TypeError('Parameter methods has to be a list.')
        if source not in [None, 1, 2]:
            raise ValueError('In Model.set_magnification_methods() ' +
                             'the parameter source, has to be 1, 2 or None.')

        if source is None:
            if isinstance(self._methods, dict):
                raise ValueError('You cannot set methods for all sources ' +
                                 'after setting them for a single source')
            self._methods = methods
        else:
            if isinstance(self._methods, list):
                raise ValueError('You cannot set methods for a single ' +
                                 'source after setting them for all sources.')
            if source > self.n_sources:
                msg = ('Cannot set methods for source {:} for model with ' +
                       'only {:} sources.')
                raise ValueError(msg.format(source, self.n_sources))
            if self._methods is None:
                self._methods = {}
            self._methods[source] = methods

    def set_default_magnification_method(self, method):
        """
        Stores information on method to be used, when no method is
        directly specified. See
        :py:class:`~MulensModel.magnificationcurve.MagnificationCurve`
        for a list of implemented methods.

        Parameters:
            method: *str*
                Name of the method to be used.

        """
        self._default_magnification_method = method

    def set_magnification_methods_parameters(self, methods_parameters):
        """
        Set additional parameters for magnification calculation methods.

        Parameters :
            methods_parameters: *dict*
                Dictionary that for method names (keys) returns dictionary
                in the form of ``**kwargs`` that are passed to given method,
                e.g., ``{'VBBL': {'accuracy': 0.005}}``.

        """
        if self.n_lenses == 1:
            methods_ok = [
                'point_source',
                'finite_source_uniform_Gould94'.lower(),
                'finite_source_uniform_Gould94_direct'.lower(),
                'finite_source_LD_Yoo04'.lower(),
                'finite_source_LD_Yoo04_direct'.lower(),
                'finite_source_uniform_Lee09'.lower(),
                'finite_source_LD_Lee09'.lower()]
        elif self.n_lenses == 2:
            methods_ok = [
                'point_source', 'quadrupole', 'hexadecapole', 'vbbl',
                'adaptive_contouring', 'point_source_point_lens']
        else:
            msg = 'wrong value of Model.n_lenses: {:}'
            raise ValueError(msg.format(self.n_lenses))

        parameters = {
            key.lower(): value for (key, value) in methods_parameters.items()}
        methods = set(parameters.keys()) - set(methods_ok)

        if len(methods):
            raise KeyError('Unknown methods provided: {:}'.format(methods))

        self._methods_parameters = parameters

    def set_limb_coeff_gamma(self, bandpass, coeff):
        """
        Store gamma limb darkening coefficient for given band. See
        also
        :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.

        Parameters :
            bandpass: *str*
                Bandpass for the coefficient you provide.

            coeff: *float*
                Value of the coefficient.

        """
        if bandpass not in self._bandpasses:
            self._bandpasses.append(bandpass)
        self._limb_darkening_coeffs.set_limb_coeff_gamma(bandpass, coeff)

    def get_limb_coeff_gamma(self, bandpass):
        """
        Get gamma limb darkening coefficient for given band.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient will be provided.

        Returns :
            gamma: *float*
                limb darkening coefficient

        """
        return self._limb_darkening_coeffs.get_limb_coeff_gamma(bandpass)

    def set_limb_coeff_u(self, bandpass, coeff):
        """
        Store u limb darkening coefficient for given band.  See also
        :py:class:`MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient you provide.

            coeff: *float*
                Value of the coefficient.

        """
        if bandpass not in self._bandpasses:
            self._bandpasses.append(bandpass)
        self._limb_darkening_coeffs.set_limb_coeff_u(bandpass, coeff)

    def get_limb_coeff_u(self, bandpass):
        """
        Get u limb darkening coefficient for given band.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient will be provided.

        Returns :
            u: *float*
                limb darkening coefficient

        """
        return self._limb_darkening_coeffs.get_limb_coeff_u(bandpass)

    def set_source_flux_ratio(self, ratio):
        """
        Sets flux ratio of sources for binary source models. If you also call
        :py:func:`set_source_flux_ratio_for_band()`, then the value set here
        will be used when: 1) no band is specified, or 2) band is specified
        but flux ratio for given band was not specified.

        Parameters :
            ratio: *float* or *None*
                The ratio of fluxes of source no. 2 to source no. 1, i.e.,
                flux_source_2/flux_source_1. Setting it to *None* removes
                the internal information, i.e., flux ratio will be fitted
                via regression (unless specific value is provided for
                bandpass).
        """
        if not isinstance(ratio, (np.float, float, type(None))):
            raise TypeError(
                'wrong type of input in Model.set_source_flux_ratio(): ' +
                'got {:}, expected float or None'.format(type(ratio)))
        if ratio is None:
            del self._source_flux_ratio_constraint[None]
        else:
            self._source_flux_ratio_constraint[None] = ratio

    def set_source_flux_ratio_for_band(self, band, ratio):
        """
        Sets flux ratio for binary source models for given band.

        Parameters :
            band: *str*
                Band for which constraint is given.

            ratio: *float*
                ratio of fluxes of source no. 2 to source no. 1, i.e.,
                flux_source_band_2/flux_source_band_1
        """
        if not isinstance(band, str):
            raise TypeError((
                'wrong type of input in ' +
                'Model.set_source_flux_ratio_for_band(): got {:}, ' +
                'expected string').format(type(band)))
        if not isinstance(ratio, (np.float, float)):
            raise TypeError((
                'wrong type of input in ' +
                'Model.set_source_flux_ratio_for_band(): got {:}, ' +
                'expected float').format(type(ratio)))
        if len(self._datasets) > 0:
            bands = [d.bandpass for d in self._datasets]
            if band not in bands:
                warnings.warn("No datasets with bandpass {:}".format(band),
                              UserWarning)

        self._source_flux_ratio_constraint[band] = ratio

    def parallax(
            self, earth_orbital=None, satellite=None, topocentric=None):
        """
        Specifies the types of the microlensing parallax that will be
        included in calculations.

        Parameters :

            earth_orbital: *boolean*, optional
                Do you want to include the effect of Earth motion about
                the Sun? Default is *False*.

            satellite: *boolean*, optional
                Do you want to include the effect of difference due to
                separation between the Earth and satellite? Note that this
                separation changes over time. Default is *False*.

            topocentric: *boolean*, optional
                Do you want to include the effect of different positions
                of observatories on the Earth? Default is *False*.
                Note that this is significant only for very high magnification
                events and if high quality datasets are analyzed.
                Hence, this effect is rarely needed. **Not Implemented yet.**

        """
        if earth_orbital is None and satellite is None and topocentric is None:
            return self._parallax
        else:
            if earth_orbital is not None:
                self._parallax['earth_orbital'] = earth_orbital
            if satellite is not None:
                self._parallax['satellite'] = satellite
            if topocentric is not None:
                self._parallax['topocentric'] = topocentric

    def get_parallax(self):
        """
        Returns *dict* that specifies the types of the microlensing parallax
        that are included in calculations.

        Returns :
            parallax: *dict*
                For keys ``'earth_orbital'``, ``'satellite'``,
                and ``'topocentric'`` returns *bool*.
        """
        return self._parallax

    def get_satellite_coords(self, times):
        """
        Get *astropy.SkyCoord* object that gives satellite positions
        for given times. see also
        :py:class:`MulensModel.satelliteskycoord.SatelliteSkyCoord`

        Parameters :
            times: *np.ndarray* or *list*
                Epochs for which satellite position is requested.

        Returns :
            satellite_skycoord: *astropy.SkyCoord*
                *SkyCoord* giving satellite positions. The parameter
                *representation* is set to 'spherical'. If
                `ephemerides_file` is not set, returns *None*.

        """
        if self.ephemerides_file is None:
            return None
        else:
            satellite_skycoords = SatelliteSkyCoord(
                 ephemerides_file=self.ephemerides_file)
            return satellite_skycoords.get_satellite_coords(times)

    def magnification(self, time, satellite_skycoord=None, gamma=0.,
                      flux_ratio_constraint=None, separate=False):
        """
        Calculate the model magnification for the given time(s).

        Parameters :
            time: *np.ndarray*, *list of floats*, or *float*
                Times for which magnification values are requested.

            satellite_skycoord: *astropy.coordinates.SkyCoord*, optional
                *SkyCoord* object that gives satellite positions. Must be
                the same length as time parameter. Use only for satellite
                parallax calculations.

            gamma: *float*, optional
                The limb darkening coefficient in gamma convention. Default is
                0 which means no limb darkening effect.

            flux_ratio_constraint:
            :py:class:`~MulensModel.mulensdata.MulensData` or *str*, optional
                Data to constrain the flux ratio of sources in binary source
                models. Can be :py:class:`~MulensModel.mulensdata.MulensData`
                instance and in that case this dataset is used to find flux
                ratio via regression and this flux ratio is applied in
                calculation of effective magnification.  If *str* is provided,
                then it indicates the bandpass and we use value set by
                :py:func:`set_source_flux_ratio_for_band()`.

            separate: *boolean*, optional
                For binary source models, return magnification of each source
                separately. Default is *False* and then only effective
                magnification is returned.

        Returns :
            magnification: *np.ndarray*
                A vector of calculated magnification values. For binary source
                models, the effective magnification is returned (unless
                *separate=True*).
        """
        allowed = (MulensData, str, type(None))
        if not isinstance(flux_ratio_constraint, allowed):
            raise TypeError(
                'The flux_ratio_constraint must be MulensData or str. If ' +
                'you want to fix it, then pass float value to ' +
                'Model.set_source_flux_ratio().\n' +
                'Got: {:}'.format(flux_ratio_constraint))
        mag = self._magnification(time, satellite_skycoord, gamma,
                                  flux_ratio_constraint, separate,
                                  same_dataset=False)
        return mag

    def _magnification(self, time, satellite_skycoord, gamma,
                       flux_ratio_constraint, separate, same_dataset):
        """
        Internal function that calculates magnification.
        """
        # Check for type
        if not isinstance(time, np.ndarray):
            if isinstance(time, (np.float, float, np.int, int)):
                time = np.array([time])
            elif isinstance(time, list):
                time = np.array(time)
            else:
                raise TypeError(
                    'time must be a float, int, list, or np.ndarray')

        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(time)

        if self.n_sources == 1:
            if flux_ratio_constraint is not None:
                raise ValueError(
                    'Model.magnification() parameter ' +
                    'flux_ratio_constraint has to be None for single source ' +
                    'models, not {:}'.format(type(flux_ratio_constraint)))
            if separate:
                raise ValueError(
                    'Model.magnification() parameter separate ' +
                    'cannot be True for single source models')
            magnification = self._magnification_1_source(
                                time, satellite_skycoord, gamma)
        elif self.n_sources == 2 and separate:
            if flux_ratio_constraint is not None:
                raise ValueError(
                    'You cannot set both flux_ratio_constraint and separate' +
                    " parameters in Model.magnification(). This doesn't make" +
                    'sense')
            magnification = self._magnification_2_sources(
                time, satellite_skycoord, gamma, flux_ratio_constraint,
                separate, same_dataset)
        elif self.n_sources == 2:
            dict_constraint = self._source_flux_ratio_constraint
            if isinstance(flux_ratio_constraint, MulensData):
                band = flux_ratio_constraint.bandpass
                if band is not None and band in dict_constraint:
                    flux_ratio_constraint = dict_constraint[band]
                elif None in dict_constraint:
                    flux_ratio_constraint = dict_constraint[None]
            elif isinstance(flux_ratio_constraint, str):
                flux_ratio_constraint = dict_constraint[flux_ratio_constraint]
            elif None in dict_constraint:
                flux_ratio_constraint = dict_constraint[None]
            magnification = self._magnification_2_sources(
                                time, satellite_skycoord, gamma,
                                flux_ratio_constraint, separate, same_dataset)
        else:
            raise ValueError('strange number of sources: {:}'.format(
                    self.n_sources))
        return magnification

    def _magnification_1_source(self, time, satellite_skycoord, gamma):
        """
        calculate model magnification for given times for model with
        a single source
        """
        magnification_curve = MagnificationCurve(
            time, parameters=self.parameters,
            parallax=self._parallax, coords=self._coords,
            satellite_skycoord=satellite_skycoord,
            gamma=gamma)
        magnification_curve.set_magnification_methods(
            self._methods, self._default_magnification_method)
        magnification_curve.set_magnification_methods_parameters(
            self._methods_parameters)

        return magnification_curve.magnification

    def _magnification_2_sources(self, time, satellite_skycoord, gamma,
                                 flux_ratio_constraint,
                                 separate, same_dataset):
        """
        calculate model magnification for given times for model with
        two sources

        flux_ratio_constraint: *float* or *MulensData*

        same_dataset: *boolean*
            If *flux_ratio_constraint* is of *MulensData* type, then is it
            the same dataset as the one for which you want magnification?
        """
        (mag_1, mag_2) = self._separate_magnifications(
                time, satellite_skycoord, gamma)

        if separate:
            return (mag_1, mag_2)

        allowed = (MulensData, float, np.float)
        if not isinstance(flux_ratio_constraint, allowed):
            raise TypeError(
                'Source flux ratio must me float or MulensData at this ' +
                'point, not {:}'.format(type(flux_ratio_constraint)))

        if isinstance(flux_ratio_constraint, (float, np.float)):
            source_flux_ratio = flux_ratio_constraint
            self._fit = None
        else:
            if same_dataset:
                if not (flux_ratio_constraint.time == time).all():
                    raise ValueError(
                        'internal problem with time vectors:\n' +
                        str(flux_ratio_constraint.time) + '\n' + str(time))
                mags = np.array([mag_1, mag_2])
            elif not same_dataset:
                (mag_1_data, mag_2_data) = self._separate_magnifications(
                    flux_ratio_constraint.time, satellite_skycoord, gamma)
                mags = np.array([mag_1_data, mag_2_data])
            else:
                raise ValueError('same_dataset must be True or False, not ' +
                                 '{:}'.format(same_dataset))
            self._fit = Fit(data=flux_ratio_constraint, magnification=mags)
            self._fit.fit_fluxes()
            f_s = self._fit.flux_of_sources(flux_ratio_constraint)
            source_flux_ratio = f_s[1] / f_s[0]

        magnification = mag_1 + mag_2 * source_flux_ratio
        magnification /= (1. + source_flux_ratio)
        return magnification

    def _separate_magnifications(self, time, satellite_skycoord, gamma):
        """
        Calculate magnification separately for each source.
        """
        kwargs = {'times': time, 'parallax': self._parallax,
                  'coords': self._coords,
                  'satellite_skycoord': satellite_skycoord, 'gamma': gamma}

        if isinstance(self._methods, dict):
            methods_1 = self._methods.get(1, None)
            methods_2 = self._methods.get(2, None)
        else:
            methods_1 = self._methods
            methods_2 = self._methods

        self._magnification_curve_1 = MagnificationCurve(
            parameters=self.parameters.source_1_parameters, **kwargs)
        self._magnification_curve_1.set_magnification_methods(
            methods_1, self._default_magnification_method)
        self._magnification_curve_1.set_magnification_methods_parameters(
            self._methods_parameters)
        mag_1 = self._magnification_curve_1.magnification

        self._magnification_curve_2 = MagnificationCurve(
            parameters=self.parameters.source_2_parameters, **kwargs)
        self._magnification_curve_2.set_magnification_methods(
            methods_2, self._default_magnification_method)
        self._magnification_curve_2.set_magnification_methods_parameters(
            self._methods_parameters)
        mag_2 = self._magnification_curve_2.magnification

        return (mag_1, mag_2)

    @property
    def data_magnification(self):
        """
        *list*

        A list of magnifications calculated for every dataset in
        :py:attr:`datasets`.
        """
        self._data_magnification = []

        fit = None
        for dataset in self._datasets:
            magnification = self.get_data_magnification(dataset)
            self._data_magnification.append(magnification)
            if (self._fit is not None and self.n_sources > 1 and
                    None not in self._source_flux_ratio_constraint):
                if fit is None:
                    fit = self._fit
                else:
                    fit.update(self._fit)
        if self.n_sources > 1:
            self._fit = fit

        return self._data_magnification

    def get_data_magnification(self, dataset):
        """
        Get the model magnification for a dataset.

        Parameters :
            dataset: :py:class:`~MulensModel.mulensdata.MulensData`
                Dataset with epochs for which magnification will be given.
                Satellite and limb darkening information is taken into
                account.

        Returns :
            magnification_vector: *np.ndarray*
                Values on magnification.

        """
        if dataset.ephemerides_file is not None:
            dataset_satellite_skycoord = dataset.satellite_skycoord
        else:
            dataset_satellite_skycoord = None

        if not self.parameters.is_finite_source():
            gamma = 0.
        else:
            if dataset.bandpass is None:
                gamma = 0.
            elif dataset.bandpass not in self.bandpasses:
                message = (
                    "Limb darkening coefficient requested for bandpass {:}, " +
                    "but not set before. We're assuming gamma=0. Use " +
                    "set_limb_coeff_gamma() or set_limb_coeff_u().")
                warnings.warn(message.format(dataset.bandpass), RuntimeWarning)
                gamma = 0.
            else:
                gamma = self._limb_darkening_coeffs.get_limb_coeff_gamma(
                    dataset.bandpass)

        if self.parameters.n_sources == 1:
            flux_ratio_constraint = None
        elif self.parameters.n_sources == 2:
            flux_ratio_constraint = dataset
        else:
            raise ValueError('Wrong number of sources')
        magnification = self._magnification(
                dataset.time, satellite_skycoord=dataset_satellite_skycoord,
                gamma=gamma, flux_ratio_constraint=flux_ratio_constraint,
                separate=False, same_dataset=True)
        return magnification

    @property
    def caustics(self):
        """
        :py:class:`~MulensModel.caustics.Caustics`

        Caustics for given model
        """
        return self._caustics

    @property
    def parameters(self):
        """
        :py:class:`~MulensModel.modelparameters.ModelParameters`

        Model parameters.
        """
        return self._parameters

    @parameters.setter
    def parameters(self, new_params):
        if isinstance(new_params, ModelParameters):
            self._parameters = new_params
        elif isinstance(new_params, dict):
            self._parameters = ModelParameters(new_params)
        else:
            raise TypeError(
                'Model.parameters must be a dictionary or ModelParameters ' +
                'object.')

    @property
    def n_lenses(self):
        """
        *int*

        number of objects in the lens system
        """
        return self._parameters.n_lenses

    @property
    def n_sources(self):
        """
        *int*

        number of luminous sources; it's possible to be 1 for xallarap model
        """
        return self._parameters.n_sources

    def is_static(self):
        """
        see :py:func:`MulensModel.modelparameters.ModelParameters.is_static()`
        """
        return self._parameters.is_static()

    @property
    def coords(self):
        """
        see :py:class:`~MulensModel.coordinates.Coordinates`
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        self._coords = Coordinates(new_value)

    @property
    def bandpasses(self):
        """
        *list*

        List of all bandpasses for which limb darkening coefficients are set.
        """
        return self._bandpasses

    @property
    def datasets(self):
        """
        *list* of :py:class:`~MulensModel.mulensdata.MulensData`

        Datasets linked to given model. Note that these can be changed by
        :py:class:`~MulensModel.event.Event` instances. This happens when
        the same model is linked to multiple
        :py:class:`~MulensModel.event.Event` instances.
        """
        if len(self._datasets) == 0:
            raise ValueError('No datasets were linked to the model')
        # The condition below should be removed, but we keep it for backward
        # compatibility, i.e., at some point if user called
        # Model.set_datasets(dataset) then Model.datasets was dataset, but for
        # Model.set_datasets([dataset]) then Model.datasets was [dataset].
        # Note that always self._datasets is of list type.
        if len(self._datasets) == 1 and self._datasets_set_single:
            return self._datasets[0]
        return self._datasets

    def set_datasets(self, datasets, data_ref=0):
        """
        Set :obj:`datasets` property

        Parameters :
            datasets: *list* of :py:class:`~MulensModel.mulensdata.MulensData`
                Datasets to be stored.

            data_ref: *int* or,
            :py:class:`~MulensModel.mulensdata.MulensData`, optional

                Reference dataset.
        """
        if isinstance(datasets, MulensData):
            self._datasets_set_single = True
            self._datasets = [datasets]
        else:
            self._datasets_set_single = False
            self._datasets = datasets
        self._data_magnification = None
        self.data_ref = data_ref

    @property
    def fit(self):
        """
        :py:class:`MulensModel.fit.Fit`

        :py:class:`MulensModel.fit.Fit` instance recently used. It gives
        access to source and blending fluxes.
        """
        return self._fit
