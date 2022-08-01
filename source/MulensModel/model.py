import warnings
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord

from MulensModel.caustics import Caustics
from MulensModel.causticswithshear import CausticsWithShear
from MulensModel.coordinates import Coordinates
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.modelparameters import ModelParameters
from MulensModel.mulensdata import MulensData
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.trajectory import Trajectory
from MulensModel.utils import Utils, PlotUtils


class Model(object):
    """
    A Model for a microlensing event with the specified parameters.

    Arguments :
        parameters: *dictionary*,
        :py:class:`~MulensModel.modelparameters.ModelParameters`

            see
            :py:class:`MulensModel.modelparameters.ModelParameters`

        :py:obj:`coords`: *str*, *astropy.SkyCoords*,
        *MulensModel.Coordinates*, optional

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

        self._limb_darkening_coeffs = LimbDarkeningCoeffs()
        self._bandpasses = []

    def __repr__(self):
        return '{0}'.format(self.parameters)

    def plot_magnification(
            self, times=None, t_range=None, t_start=None, t_stop=None, dt=None,
            n_epochs=None, subtract_2450000=False, subtract_2460000=False,
            satellite_skycoord=None, gamma=None, source_flux_ratio=None,
            flux_ratio_constraint=None,
            **kwargs):
        """
        Plot the model magnification curve.

        Keywords :
            see :py:func:`plot_lc()`

            gamma:
                see :py:func:`get_magnification()`

            satellite_skycoord:
                see :py:func:`plot_trajectory()`

            source_flux_ratio: *float*
                If the model has two sources, source_flux_ratio is the ratio of
                source_flux_2 / source_flux_1

            flux_ratio_constraint: DEPRECATED. Use source_flux_ratio instead.

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.

        """
        if flux_ratio_constraint is not None:
            warnings.warn(
                'flux_ratio_constraint will be deprecated. Use ' +
                'source_flux_ratio instead')
            source_flux_ratio = flux_ratio_constraint

        if self.n_sources > 1 and source_flux_ratio is None:
            raise ValueError(
                'For binary source model you have to provide ' +
                'source_flux_ratio. Note that plotted magnification will ' +
                'be the effective magnification of the two sources.')

        if 'fit_blending' in kwargs:
            raise AttributeError(
                'fit_blending is deprecated. See Event() class instead.')

        self._check_gamma_for_2_sources(gamma)

        if times is None:
            times = self.set_times(
                t_range=t_range, t_start=t_start, t_stop=t_stop, dt=dt,
                n_epochs=n_epochs)

        subtract = PlotUtils.find_subtract(subtract_2450000=subtract_2450000,
                                           subtract_2460000=subtract_2460000)

        if satellite_skycoord is not None:
            if isinstance(satellite_skycoord, SatelliteSkyCoord):
                satellite = satellite_skycoord.get_satellite_coords(times)
            elif not isinstance(satellite_skycoord, SkyCoord):
                raise TypeError('Wrong type of satellite_skycoord in ' +
                                'Model.plot_magnification()')
        else:
            satellite = None

        magnification = self.get_magnification(
            times, satellite_skycoord=satellite, gamma=gamma,
            source_flux_ratio=source_flux_ratio)

        self._plt_plot(times-subtract, magnification, kwargs)
        plt.ylabel('Magnification')
        plt.xlabel(
            PlotUtils.find_subtract_xlabel(
                subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000))

    def get_lc(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, source_flux=None, blend_flux=None,
            source_flux_ratio=None, gamma=None, bandpass=None):
        """
        Calculate model light curve in magnitudes.

        Keywords :
            times: [*float*, *list*, *numpy.ndarray*]
                a list of times at which to plot the magnifications

            t_range, t_start, t_stop, dt, n_epochs: see :py:func:`set_times`

            source_flux: *float* or *list*
                Explicitly specify the source flux(es) in a
                system where flux = 1 corresponds to
                :obj:`MulensModel.utils.MAG_ZEROPOINT` (= 22 mag). If the model
                has n_source > 1, source_flux may be specified as a list: one
                value for each source. Alternatively, if source_flux is
                specified as a float, source_flux_ratio should also be
                specified. Then, source_flux is taken to be the flux of the
                first source, and the other source fluxes are derived using
                source_flux_ratio.

            blend_flux: *float*
                Explicitly specify the blend flux in a
                system where flux = 1 corresponds to
                :obj:`MulensModel.utils.MAG_ZEROPOINT` (= 22 mag).

            source_flux_ratio: *float*, Optional
                If the model has two sources, source_flux_ratio is the ratio of
                source_flux_2 / source_flux_1.

            gamma, bandpass:
                see :py:func:`get_magnification()`

        Returns :
            magnification: *numpy.ndarray*
                Magnification values for each epoch.
        """
        fluxes = self._parse_fluxes_for_get_lc(
            source_flux, source_flux_ratio, blend_flux)
        (source_flux, source_flux_ratio, blend_flux) = fluxes

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)
        self._check_gamma_for_2_sources(gamma)

        magnitudes = self._get_lc(
            times=times, t_range=t_range, t_start=t_start, t_stop=t_stop,
            dt=dt, n_epochs=n_epochs, gamma=gamma, source_flux=source_flux,
            blend_flux=blend_flux, return_times=False)

        return magnitudes

    def _check_gamma_for_2_sources(self, gamma):
        """
        Check if the user tries to use limb darkening for binary source model
        with finite source effect of both sources. If that is the case,
        then raise exception.
        The gamma value of *None* or *0* indicates that limb-darkening is off.
        """
        if gamma is None or float(gamma) == 0.0:
            return
        if self.n_sources == 1:
            return

        is_finite_1 = self._parameters.source_1_parameters.is_finite_source()
        is_finite_2 = self._parameters.source_2_parameters.is_finite_source()
        if is_finite_1 and is_finite_2:
            raise NotImplementedError(
                "You're requesting binary source model with both sources " +
                "showing finite source effect and you're specifying " +
                "limb-darkening coefficient. This is not yet implemented.")

    def _parse_fluxes_for_get_lc(self, source_flux, source_flux_ratio,
                                 blend_flux):
        """
        Parsing of fluxes to be used in get_lc/plot_lc.
        """
        if source_flux is None:
            raise ValueError("You must provide a value for source_flux.")
        elif (isinstance(source_flux, float) and self.n_sources > 1):
            if source_flux_ratio is None:
                raise ValueError(
                    "Either source_flux should be a list or " +
                    "source_flux_ratio should be specified.\n" +
                    "source_flux = {0}\n".format(source_flux) +
                    "n_sources = {0}\n".format(self.n_sources) +
                    "source_flux_ratio = {0}".format(source_flux_ratio))
            else:
                source_flux = [source_flux]
                source_flux.append(source_flux[0] * source_flux_ratio)

        if blend_flux is None:
            warnings.warn(
                'No blend_flux not specified. Assuming blend_flux = zero.')
            blend_flux = 0.

        return (source_flux, source_flux_ratio, blend_flux)

    def _get_lc(self, times, t_range, t_start, t_stop, dt, n_epochs, gamma,
                source_flux, blend_flux, return_times=False):
        """
        calculate magnitudes without making checks on input parameters

        source_flux is a *float* (for single source model) or
        an iterable (for multiple sources)
        """
        if times is None:
            times = self.set_times(
                t_range=t_range, t_start=t_start, t_stop=t_stop, dt=dt,
                n_epochs=n_epochs)
        elif isinstance(times, list):
            times = np.array(times)

        if self.n_sources == 1:
            magnification = self.get_magnification(times, gamma=gamma)
            flux = source_flux * magnification + blend_flux
        else:
            magnification = self.get_magnification(times, separate=True)
            flux = None
            for i in range(self.n_sources):
                if flux is None:
                    flux = source_flux[i] * magnification[i]
                else:
                    flux += source_flux[i] * magnification[i]

            flux += blend_flux

        magnitudes = Utils.get_mag_from_flux(flux)

        if return_times:
            return (times, magnitudes)
        else:
            return magnitudes

    def plot_lc(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, source_flux=None, blend_flux=None,
            source_flux_ratio=None, gamma=None, bandpass=None,
            subtract_2450000=False, subtract_2460000=False,
            data_ref=None, flux_ratio_constraint=None,
            fit_blending=None, f_source=None, f_blend=None,
            **kwargs):
        """
        Plot the model light curve in magnitudes.

        Keywords :
            times: [*float*, *list*, *numpy.ndarray*]
                a list of times at which to plot the magnifications

            t_range, t_start, t_stop, dt, n_epochs:
                see :py:func:`set_times()`

            source_flux, blend_flux, source_flux_ratio:
                see :py:func:`get_lc()`

            gamma, bandpass:
                see :py:func:`get_magnification()`

            subtract_2450000, subtract_2460000: *boolean*, optional
                If True, subtracts 2450000 or 2460000 from the time
                axis to get more human-scale numbers. If using, make
                sure to also set the same settings for all other
                plotting calls (e.g. :py:func:`plot_data()`)

            data_ref: DEPRECATED
                Specify source_flux and blend_flux instead or use plotting
                functions in py:class:`~MulensModel.Event()`

            flux_ratio_constraint: DEPRECATED
                Use source_flux_ratio instead

            fit_blending: DEPRECATED
                Use py:class:`~MulensModel.Event()` for fitting.

            f_source, f_blend: DEPRECATED
                use *source_flux* or *blend_flux* instead.

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.
        """

        if flux_ratio_constraint is not None:
            warnings.warn(
                'flux_ratio_constraint will be deprecated. Use ' +
                'source_flux_ratio instead')
            source_flux_ratio = flux_ratio_constraint

        if data_ref is not None:
            raise AttributeError(
                'data_ref keyword has been deprecated. Specify source_flux ' +
                'and blend_flux instead or use plotting functions in Event().')

        if fit_blending is not None:
            raise AttributeError(
                'fit_blending keyword has been deprecated. Use Event() ' +
                'instead.')

        if f_source is not None:
            warnings.warn(
                'f_source will be deprecated. Use source_flux instead')
            source_flux = f_source

        if f_blend is not None:
            warnings.warn(
                'f_blend will be deprecated. Use blend_flux instead')
            blend_flux = f_blend

        fluxes = self._parse_fluxes_for_get_lc(
            source_flux, source_flux_ratio, blend_flux)
        (source_flux, source_flux_ratio, blend_flux) = fluxes

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)
        self._check_gamma_for_2_sources(gamma)

        (times, magnitudes) = self._get_lc(
            times=times, t_range=t_range, t_start=t_start, t_stop=t_stop,
            dt=dt, n_epochs=n_epochs, gamma=gamma, source_flux=source_flux,
            blend_flux=blend_flux, return_times=True)

        subtract = PlotUtils.find_subtract(subtract_2450000=subtract_2450000,
                                           subtract_2460000=subtract_2460000)

        self._plt_plot(times-subtract, magnitudes, kwargs)
        plt.ylabel('Magnitude')
        plt.xlabel(
            PlotUtils.find_subtract_xlabel(
                subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000))

        (ymin, ymax) = plt.gca().get_ylim()
        if ymax > ymin:
            plt.gca().invert_yaxis()

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

        # check if covergence_K and shear_G are in parameters
        if self.parameters.is_external_mass_sheet:
            self._caustics = CausticsWithShear(
                q=self.parameters.q, s=s,
                convergence_K=self.parameters.convergence_K,
                shear_G=self.parameters.shear_G)
        else:
            self._caustics = Caustics(q=self.parameters.q, s=s)

    def plot_trajectory(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, caustics=False,
            arrow=True, satellite_skycoord=None, arrow_kwargs=None,
            show_data=None, **kwargs):
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

            show_data: DEPRECATED
                Use py:class:`~MulensModel.Event()` for plotting data with
                models.

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
        if show_data is not None:
            raise AttributeError(
                'show_data is deprecated. datasets are no longer part of ' +
                'Model. See Event.plot_source_for_datasets() instead.')

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
        times = np.atleast_1d(times)

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

    def get_trajectory(self, times):
        """
        Get the source trajectory for the given set of times.

        Parameters :
            times:  *np.ndarray*, *list of floats*, or *float*
                Times for which magnification values are requested.

        Returns : A `:py:class:`~MulensModel.trajectory.Trajectory` object.

        """
        kwargs_ = {
            'times': times, 'parallax': self._parallax, 'coords': self._coords,
            'satellite_skycoord': self.get_satellite_coords(times)}
        return Trajectory(parameters=self.parameters, **kwargs_)

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

    def _get_limb_coeff_gamma(self, bandpass, gamma):
        """
        Get gamma from either bandpass or gamma
        """
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

        return gamma

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

    def get_magnification(self, time, satellite_skycoord=None, gamma=None,
                          bandpass=None, source_flux_ratio=None,
                          separate=False, flux_ratio_constraint=None):
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
                The limb-darkening coefficient in gamma convention. Default is
                0 which means no limb darkening effect.

            bandpass: *str*, optional
                The bandpass for setting the limb-darkening coefficient.
                Expects that you have used :py:func:`set_limb_coeff_gamma()` or
                :py:func:`set_limb_coeff_u()`. Only ONE of 'gamma' or
                'bandpass' may be specified.

            source_flux_ratio: *float*
                If the model has two sources, source_flux_ratio is the ratio of
                source_flux_2 / source_flux_1

            separate: *boolean*, optional
                For binary source models, return magnification of each source
                separately. Default is *False* and then only effective
                magnification is returned.

            flux_ratio_constraint: DEPRECATED
                Use source_flux_ratio instead.

        Returns :
            magnification: *np.ndarray*
                A vector of calculated magnification values. For binary source
                models, the effective magnification is returned (unless
                *separate=True*).
        """
        if flux_ratio_constraint is not None:
            warnings.warn(
                'flux_ratio_constraint will be deprecated. Use ' +
                'source_flux_ratio instead.')
            if isinstance(flux_ratio_constraint, float):
                source_flux_ratio = flux_ratio_constraint
            elif isinstance(flux_ratio_constraint, MulensData):
                raise AttributeError(
                    'The ability to set flux_ratio_constraint with a dataset' +
                    'is deprecated. Use a float with source_flux_ratio ' +
                    'instead.')
            else:
                raise ValueError(
                    'Wrong type for flux_ratio_constraint. Use a float with '
                    'source_flux_ratio instead.')
        elif source_flux_ratio is not None:
            if not isinstance(source_flux_ratio, float):
                raise TypeError(
                    'source_flux_ratio should be a float. Got: {:}'.format(
                        source_flux_ratio))

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)
        self._check_gamma_for_2_sources(gamma)

        if self.n_sources > 1:
            if (source_flux_ratio is None) and (separate is False):
                raise ValueError(
                    'For 2 sources either source_flux_ratio should be set or' +
                    ' separate=True. \n' +
                    'separate: {0}\n'.format(separate) +
                    'source_flux_ratio: {0}'.format(source_flux_ratio))

        magnification = self._get_magnification(
            time, satellite_skycoord, gamma, source_flux_ratio, separate)

        return magnification

    def magnification(self, *args, **kwargs):
        """
        DEPRECATED

        Use :py:func:`get_magnification()` instead.
        """
        warnings.warn('magnification() will be deprecated in ' +
                      'favor of get_magnification()')
        return self.get_magnification(*args, **kwargs)

    def _get_magnification(self, time, satellite_skycoord, gamma,
                           source_flux_ratio, separate):
        """
        Internal function that calculates magnification.
        """
        time = np.atleast_1d(time)

        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(time)

        if self.n_sources == 1:
            if source_flux_ratio is not None:
                raise ValueError(
                    'Model.get_magnification() parameter ' +
                    'flux_ratio_constraint has to be None for single source ' +
                    'models, not {:}'.format(source_flux_ratio))
            elif separate:
                raise ValueError(
                    'Model.get_magnification() parameter separate ' +
                    'cannot be True for single source models')
            else:
                magnification = self._magnification_1_source(
                    time, satellite_skycoord, gamma)

        elif self.n_sources == 2:
            magnification = self._magnification_2_sources(
                time, satellite_skycoord, gamma, source_flux_ratio,
                separate)
        else:
            raise ValueError(
                'Only 1 or 2 sources is implemented. Number of sources: ' +
                '{:}'.format(self.n_sources))

        if np.sum(np.isnan(magnification)) > 0:
            fmt = ("EPOCHS:\n{:}\nMODEL:\n{:}Something went wrong with " +
                   "calculating of magnification for the above model. " +
                   "For all above epochs magnifications are NaN.")
            msg = fmt.format(time[np.isnan(magnification)], self.__repr__())
            raise ValueError(msg)
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

        return magnification_curve.get_magnification()

    def _magnification_2_sources(
            self, time, satellite_skycoord, gamma, source_flux_ratio,
            separate):
        """
        calculate model magnification for given times for model with
        two sources

        source_flux_ratio: *float*
        separate: *bool*
        """
        if separate and (source_flux_ratio is not None):
            raise ValueError(
                'You cannot set both source_flux_ratio and separate' +
                " parameters in Model.get_magnification(). This doesn't " +
                'make sense')

        (mag_1, mag_2) = self._separate_magnifications(
            time, satellite_skycoord, gamma)

        if separate:
            return (mag_1, mag_2)
        else:
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
        mag_1 = self._magnification_curve_1.get_magnification()

        self._magnification_curve_2 = MagnificationCurve(
            parameters=self.parameters.source_2_parameters, **kwargs)
        self._magnification_curve_2.set_magnification_methods(
            methods_2, self._default_magnification_method)
        self._magnification_curve_2.set_magnification_methods_parameters(
            self._methods_parameters)
        mag_2 = self._magnification_curve_2.get_magnification()

        return (mag_1, mag_2)

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

# ---- DEPRECATED PROPERTIES AND FUNCTIONS---- #
    def reset_plot_properties(self):
        """
        DEPRECATED

        Resets internal plotting properties of all attached datasets.
        """
        raise AttributeError(
            'reset_plot_properties is deprecated. datasets are ' +
            'no longer part of Model(). Use Event() instead.')

    def plot_data(
            self, data_ref=None, show_errorbars=None, show_bad=None,
            color_list=None, marker_list=None, size_list=None,
            label_list=None, alpha_list=None, zorder_list=None,
            subtract_2450000=False, subtract_2460000=False, **kwargs):
        """
        DEPRECATED

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
        raise AttributeError(
            'plot_data is deprecated. datasets are no ' +
            'longer part of Model(). Use Event() instead.')

    def plot_residuals(
            self, show_errorbars=None,
            color_list=None, marker_list=None, size_list=None,
            label_list=None, alpha_list=None, zorder_list=None,
            data_ref=None, subtract_2450000=False, subtract_2460000=False,
            show_bad=None, **kwargs):
        """
        DEPRECATED

        Plot the residuals (in magnitudes) of the model.

        For explanation of keywords, see doctrings in
        :py:func:`plot_data()`. Note different order of keywords.
        """

        raise AttributeError(
            'plot_residuals is deprecated. datasets are no ' +
            'longer part of Model(). Use Event().fits[data_ref]' +
            '.get_residuals() instead.')

    def get_residuals(self, data_ref=None, type='mag', data=None):
        """
        DEPRECATED

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
        raise AttributeError(
            'get_residuals is deprecated. datasets are no ' +
            'longer part of Model(). Use Event() instead.')

    def plot_source_for_datasets(self, **kwargs):
        """
        DEPRECATED

        Plot source positions for all linked datasets. Colors used for
        each dataset are the same used for plotting photometry.

        Parameters:
            ``**kwargs``:
                see :py:func:`plot_source`
        """
        raise AttributeError(
            'plot_source_for_datasets is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead.')

    def get_ref_fluxes(self, data_ref=None, fit_blending=None):
        """
        DEPRECATED

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
        raise AttributeError(
            'get_ref_fluxes is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead.')

    def set_source_flux_ratio(self, ratio):
        """
        DEPRECATED

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
        raise AttributeError(
            'set_source_flux_ratio is deprecated. Fluxes are not ' +
            'intrinsic to Model(). Set them explicitly when ' +
            'needed or use Event().')

    def set_source_flux_ratio_for_band(self, band, ratio):
        """
        DEPRECATED

        Sets flux ratio for binary source models for given band.

        Parameters :
            band: *str*
                Band for which constraint is given.

            ratio: *float*
                ratio of fluxes of source no. 2 to source no. 1, i.e.,
                flux_source_band_2/flux_source_band_1
        """
        raise AttributeError(
            'set_source_flux_ratio_for_band is deprecated. Fluxes are not ' +
            'intrinsic to Model(). Set them explicitly when ' +
            'needed or use Event().')

    @property
    def data_magnification(self):
        """
        DEPRECATED

        *list*

        A list of magnifications calculated for every dataset in
        :py:attr:`datasets`.
        """
        raise AttributeError(
            'data_magnification is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead.')

    def get_data_magnification(self, dataset):
        """
        DEPRECATED

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
        raise AttributeError(
            'get_data_magnification is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead or use' +
            'Model.get_magnification() for the relevant times.')

    @property
    def datasets(self):
        """
        DEPRECATED

        *list* of :py:class:`~MulensModel.mulensdata.MulensData`

        Datasets linked to given model. Note that these can be changed by
        :py:class:`~MulensModel.event.Event` instances. This happens when
        the same model is linked to multiple
        :py:class:`~MulensModel.event.Event` instances.
        """
        raise AttributeError(
            'datasets is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead.')

    def set_datasets(self, datasets, data_ref=0):
        """
        DEPRECATED

        Set :obj:`datasets` property

        Parameters :
            datasets: *list* of :py:class:`~MulensModel.mulensdata.MulensData`
                Datasets to be stored.

            data_ref: *int* or,
            :py:class:`~MulensModel.mulensdata.MulensData`, optional

                Reference dataset.
        """
        raise AttributeError(
            'set_datasets is deprecated. datasets ' +
            'are no longer part of Model(). Use Event() instead.')
