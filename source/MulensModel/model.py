import warnings
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord

from MulensModel.causticsbinary import CausticsBinary
from MulensModel.causticspointwithshear import CausticsPointWithShear
from MulensModel.causticsbinarywithshear import CausticsBinaryWithShear
from MulensModel.coordinates import Coordinates
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.modelparameters import ModelParameters
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

    This will provide information on parameter values, coordinates,
    methods used for magnification calculations, and
    limb-darkening coefficients.
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

        self._limb_darkening_coeffs = [LimbDarkeningCoeffs() for _ in range(self.n_sources)]
        self._bandpasses = []

    def __repr__(self):
        out = '{0}'.format(self.parameters)
        if self.coords is not None:
            out += '\ncoords: {0}'.format(self.coords)

        out += '\ndefault magnification method: {0}'.format(
            self._default_magnification_method)
        if self._methods is not None:
            out += '\nother magnification methods: {0}'.format(self._methods)

        if len(self.bandpasses) > 0:
            out += '\nlimb-darkening coeffs (gamma): {0}'.format(
                self._limb_darkening_coeffs)

        return out

    def __getattr__(self, item):
        return object.__getattribute__(self, item)

    def plot_magnification(
            self, times=None, t_range=None, t_start=None, t_stop=None, dt=None,
            n_epochs=None, subtract_2450000=False, subtract_2460000=False,
            satellite_skycoord=None, gamma=None, source_flux_ratio=None,
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

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.

        """
        if self.n_sources > 1 and source_flux_ratio is None:
            raise ValueError(
                'For multi-source model you have to provide ' +
                'source_flux_ratio. Note that plotted magnification will ' +
                'be the effective magnification of the sources.')

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
            magnitudes: *numpy.ndarray*
                Magnitude values for each epoch.
        """
        fluxes = self._parse_fluxes_for_get_lc(
            source_flux, source_flux_ratio, blend_flux)
        (source_flux, source_flux_ratio, blend_flux) = fluxes

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)

        magnitudes = self._get_lc(
            times=times, t_range=t_range, t_start=t_start, t_stop=t_stop,
            dt=dt, n_epochs=n_epochs, gamma=gamma, source_flux=source_flux,
            blend_flux=blend_flux, return_times=False)

        return magnitudes

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
                if isinstance(source_flux_ratio, (float)):
                    source_flux_ratio = [source_flux_ratio]

                source_flux = [source_flux]
                for i in range(0, self.n_sources-1):
                    source_flux.append(source_flux[0] * source_flux_ratio[i])

        if blend_flux is None:
            warnings.warn(
                'No blend_flux not specified. Assuming blend_flux = zero.')
            blend_flux = 0.

        return (source_flux, source_flux_ratio, blend_flux)

    def _get_lc(self, times, t_range, t_start, t_stop, dt, n_epochs, gamma,
                source_flux, blend_flux, return_times=False, phot_fmt="mag"):
        """
        calculate magnitude or flux without making checks on input parameters

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

        return self._return_mag_or_flux(times, flux, return_times, phot_fmt)

    def _return_mag_or_flux(self, times, flux, return_times, phot_fmt):
        """
        Obtain what is returned in function _get_lc, where phot_fmt and
        return_times are explicitly given
        """
        if phot_fmt == 'mag':
            mag_or_flux = Utils.get_mag_from_flux(flux)
        elif phot_fmt == 'flux':
            mag_or_flux = flux
        else:
            raise ValueError(
                'phot_fmt must be one of "mag", "flux", or "scaled". Your ' +
                'value: {0}'.format(phot_fmt))

        if return_times:
            return (times, mag_or_flux)
        else:
            return mag_or_flux

    def plot_lc(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, source_flux=None, blend_flux=None,
            source_flux_ratio=None, gamma=None, bandpass=None,
            subtract_2450000=False, subtract_2460000=False, phot_fmt="mag",
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

            phot_fmt: *str*
                Specifies whether the photometry is plotted in magnitude or
                flux space. Accepts either 'mag' or 'flux'. Default = 'mag'.

            ``**kwargs``:
                any arguments accepted by :py:func:`matplotlib.pyplot.plot()`.
        """

        fluxes = self._parse_fluxes_for_get_lc(
            source_flux, source_flux_ratio, blend_flux)
        (source_flux, source_flux_ratio, blend_flux) = fluxes

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)

        (times, mag_or_flux) = self._get_lc(
            times=times, t_range=t_range, t_start=t_start, t_stop=t_stop,
            dt=dt, n_epochs=n_epochs, gamma=gamma, source_flux=source_flux,
            blend_flux=blend_flux, return_times=True, phot_fmt=phot_fmt)

        subtract = PlotUtils.find_subtract(subtract_2450000=subtract_2450000,
                                           subtract_2460000=subtract_2460000)

        self._plt_plot(times-subtract, mag_or_flux, kwargs)
        self._plot_axes(phot_fmt, subtract_2450000, subtract_2460000)

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

    def _plot_axes(self, phot_fmt, subtract_2450000, subtract_2460000):
        """
        Adjust axes labels and ranges, given the inputs phot_fmt and subtract
        """
        if phot_fmt == 'mag':
            plt.ylabel('Magnitude')
        elif phot_fmt == 'flux':
            plt.ylabel('Flux')
        plt.xlabel(
            PlotUtils.find_subtract_xlabel(
                subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000))

        (ymin, ymax) = plt.gca().get_ylim()
        if ymax > ymin and phot_fmt == 'mag':
            plt.gca().invert_yaxis()

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
        mass_sheet = self.parameters.is_external_mass_sheet_with_shear
        if self.n_lenses == 1 and not mass_sheet:
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
        if self.n_lenses == 1:
            self._update_caustics_single_lens()
        elif self.n_lenses == 2:
            self._update_caustics_binary_lens(epoch)
        else:
            raise ValueError('updating triple lens caustics not yet coded')

    def _update_caustics_single_lens(self):
        """
        Update self._caustics for single lens.
        """
        convergence_K = self.parameters.parameters.get('convergence_K', 0)
        shear_G = self.parameters.parameters.get('shear_G', complex(0, 0))

        self._caustics = CausticsPointWithShear(
            convergence_K=convergence_K, shear_G=shear_G)

    def _update_caustics_binary_lens(self, epoch):
        """
        Update self._caustics for binary lens.
        """
        if epoch is None:
            s = self.parameters.s
        else:
            s = self.parameters.get_s(epoch)

        if self._caustics is not None:
            if s == self._caustics.s and self.parameters.q == self._caustics.q:
                return

        if not self.parameters.is_external_mass_sheet:
            self._caustics = CausticsBinary(q=self.parameters.q, s=s)
        else:
            convergence_K = self.parameters.parameters.get('convergence_K', 0.)
            shear_G = self.parameters.parameters.get('shear_G', complex(0, 0))
            self._caustics = CausticsBinaryWithShear(
                q=self.parameters.q, s=s, shear_G=shear_G,
                convergence_K=convergence_K)

    def plot_trajectory(
            self, times=None, t_range=None, t_start=None, t_stop=None,
            dt=None, n_epochs=None, caustics=False,
            arrow=True, satellite_skycoord=None, arrow_kwargs=None,
            **kwargs):
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
        elif self.n_sources >= 2:
            for i in range(self.n_sources):
                self._plot_single_trajectory(
                    times, self.parameters.__getattr__('source_{0}_parameters'.format(i+1)),
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
            self._plot_arrow(times, trajectory, kwargs, arrow_kwargs)

    def _plot_arrow(self, times, trajectory, kwargs, arrow_kwargs):
        """
        Plot arrow for given trajectory.
        """
        width_scaling_factor = 0.01
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
        width = width**.5 * width_scaling_factor

        color = kwargs.get('color', 'black')
        kwargs_ = {'width': width, 'color': color, 'lw': 0, 'zorder': -np.inf}
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
        elif self.n_sources >= 2:
            for i in range(self.n_sources):
                trajectory = Trajectory(
                    parameters=self.parameters.__getattr__('source_{0}_parameters'.format(i+1)), **kwargs_)
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

    def get_trajectory(self, times, satellite_skycoord=None):
        """
        Get the trajectory of the source for the given set of times.
        For multi-source models, one trajectory per source is returned.

        Parameters :
            times:  *np.ndarray*, *list of floats*, or *float*
                Epochs for which source positions are requested.

            satellite_skycoord: *astropy.SkyCoord*
                Allows the user to specify that the trajectory is calculated
                for a satellite. If *astropy.SkyCoord* object is provided,
                then these are satellite positions for all epochs.
                See also :py:func:`get_satellite_coords()`

        Returns :
            trajectories: `:py:class:`~MulensModel.trajectory.Trajectory` object or a *list* of them
                Single object for single source model, a *list* otherwise.
        """
        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(times)

        kwargs_ = {
            'times': times, 'parallax': self._parallax, 'coords': self._coords,
            'satellite_skycoord': satellite_skycoord}
        if self.n_sources == 1:
            return Trajectory(parameters=self.parameters, **kwargs_)
        elif self.n_sources == 2:
            trajectories = []
            for i in range(self.n_sources):
                trajectory = Trajectory(
                    parameters=self.parameters.__getattr__('source_{0}_parameters'.format(i + 1)), **kwargs_)
                trajectories.append(trajectory)

            return trajectories
        else:
            raise NotImplementedError(
                "only 1 or 2 sources allowed here at this point")

    def set_times(
            self, t_range=None, t_start=None, t_stop=None, dt=None,
            n_epochs=1000):
        """
        Return a list of times. If no keywords are specified, default
        is 1000 epochs from [:math:`t_0 - 1.5 * t_E`, :math:`t_0 + 1.5 * t_E`]
        range.

        For multi-source models, respectively, minimum and maximum of
        `t_0_N` values are used.

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
                t_0 = np.min(
                    [self.parameters.__getattr__('source_{0}_parameters'.format(i+1)).t_0
                     for i in range(self.n_sources)])

            t_start = t_0 - (n_tE * self.parameters.t_E)

        if t_stop is None:
            if self.n_sources == 1:
                t_0 = self.parameters.t_0
            else:
                t_0 = np.max(
                    [self.parameters.__getattr__('source_{0}_parameters'.format(i+1)).t_0
                     for i in range(self.n_sources)])

            t_stop = t_0 + (n_tE * self.parameters.t_E)

        if dt is None:
            if n_epochs is None:
                n_epochs = 1000

            n_epochs -= 1
            dt = (t_stop - t_start) / float(n_epochs)

        out = np.arange(t_start, t_stop + dt, dt)
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

            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts 1, 2, or
                *None* (i.e., all sources). Default is *None*
        """
        if not isinstance(methods, list):
            raise TypeError('Parameter methods has to be a list.')
        if (not isinstance(source, (int))) and (source is not None):
            raise ValueError("In Model.set_magnification_methods() the parameter 'source' has to be *int* or None.")

        self._check_methods(methods, source)

        if (source is None) or (self.n_sources == 1):
            if isinstance(self._methods, dict):
                raise ValueError('You cannot set methods for all sources after setting them for a single source')

            self._methods = methods
        else:
            if isinstance(self._methods, list):
                raise ValueError('You cannot set methods for a single source after setting them for all sources.')

            if source > self.n_sources:
                msg = ('Cannot set methods for source {:} for model with only {:} sources.')
                raise ValueError(msg.format(source, self.n_sources))

            if self._methods is None:
                self._methods = {}

            self._methods[source] = methods

    def _check_methods(self, methods, source):
        """
        Check consistency of methods:
        - are finite source methods used for point sources?
        """
        used_methods = set(methods[1::2])
        allowed = set(['point_source', 'point_source_point_lens'])
        difference = used_methods - allowed
        if len(difference) == 0:
            return

        fmt = ('It is impossible to use finite source method for '
               'a point source {:}: {:}')
        if self.n_sources == 1:
            if not self.parameters.is_finite_source():
                raise ValueError(fmt.format("", difference))
        elif self.n_sources >= 2:
            for i in range(self.n_sources):
                if source in [i+1, None]:
                    if not self.parameters.__getattr__('source_{0}_parameters'.format(i+1)).is_finite_source():
                        raise ValueError(fmt.format("no. {0}".format(i+1), difference))

    def _check_limb_darkening(self, methods):
        """
        Check if limb darkening is used without methods that allow it.

        Parameters :
            methods: *list*
                List of methods used for magnification calculation.
        """
        forbidden = {
            "finite_source_uniform_Gould94",
            "finite_source_uniform_Gould94_direct",
            "finite_source_uniform_WittMao94",
            "finite_source_uniform_Lee09",
        }

        methods = list(methods.values()) if isinstance(methods, dict) else [methods]
        if methods == [None]:
            methods = [[self.default_magnification_method]]

        for (i, method) in enumerate(methods):
            method = {m for m in method if isinstance(m, str)}
            method.add(self.default_magnification_method)
            for bandpass in self._bandpasses:
                gamma = self.get_limb_coeff_gamma(bandpass, source=i+1)
                if gamma != 0 and not (method - forbidden):
                    raise ValueError("Limb darkening requires at least one method that includes LD.")

    def get_magnification_methods(self, source=None):
        """
        Gets methods used for magnification calculation. See
        :py:func:`set_magnification_methods`

        Parameters :
            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts 1, 2, or
                *None* (i.e., all sources). Default is *None*.

        Returns :
            methods *list* or *dict*
                Methods used in a form of a *list* or *dict* of *lists*.
        """
        if (source is None):
            return self.methods
        elif (self.n_sources == 1):
            if source > 1:
                raise IndexError(
                    'Your model only has 1 source, but you requested ' +
                    'magnification methods for source {0}'.format(source))
            else:
                return self.methods

        else:
            return self.methods[source]

    @property
    def methods(self):
        """
        *list*

        List of methods used for magnification calculation (single source) or
        *dict* of such lists (multiple sources).
        """
        return self._methods

    @property
    def default_magnification_method(self):
        """
        Stores information on method to be used, when no method is
        directly specified. See
        :py:class:`~MulensModel.magnificationcurve.MagnificationCurve`
        for a list of implemented methods.

        Parameters:
            method: *str*
                Name of the method to be used.

        """
        return self._default_magnification_method

    @default_magnification_method.setter
    def default_magnification_method(self, new_method):
        self._default_magnification_method = new_method

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
            methods_all_str = (
                'point_source finite_source_uniform_Gould94 '
                'finite_source_uniform_Gould94_direct '
                'finite_source_uniform_WittMao94 finite_source_LD_WittMao94 '
                'finite_source_LD_Yoo04 finite_source_LD_Yoo04_direct '
                'finite_source_uniform_Lee09 finite_source_LD_Lee09')
        elif self.n_lenses == 2:
            methods_all_str = ('point_source quadrupole hexadecapole vbbl '
                               'adaptive_contouring point_source_point_lens')
        else:
            msg = 'wrong value of Model.n_lenses: {:}'
            raise ValueError(msg.format(self.n_lenses))

        parameters = {
            key.lower(): value for (key, value) in methods_parameters.items()}
        methods_all = set([m.lower() for m in methods_all_str.split()])
        methods = set(parameters.keys()) - methods_all

        if len(methods):
            raise KeyError('Unknown methods provided: {:}'.format(methods))
        self._check_magnification_methods_parameters(methods_parameters)

        self._methods_parameters = parameters

    def _check_magnification_methods_parameters(self, methods_parameters):
        """
        Check if the provided kwargs are valid for the given method.
        """
        msg = "{:} method allows {:} parameters, but got '{:}'."
        allowed = {'vbbl': ['accuracy'],
                   'adaptive_contouring': ['accuracy', 'ld_accuracy']}

        for method, kwargs in methods_parameters.items():
            method = method.lower()
            if method in allowed:
                invalid = set(kwargs) - set(allowed[method])
                if invalid:
                    raise KeyError(msg.format(method, allowed[method], invalid))

    def get_magnification_methods_parameters(self, method):
        """
        Get additional parameters for a specific magnification calculation
        method or methods.

        Parameters :
            method: *str*, *list*
                Name of method or a list of the names for which parameters
                will be returned.

        Returns :
            method_parameters: *dict*
                see :py:func:`set_magnification_methods_parameters`
        """
        if isinstance(method, (str)):
            parameters = {
                method.lower(): self._methods_parameters[method.lower()]}
        else:
            parameters = {key.lower(): self._methods_parameters[key.lower()]
                          for key in method}

        return parameters

    def set_limb_coeff_gamma(self, bandpass, coeff, source=None):
        """
        Store gamma limb darkening coefficient for given band. See
        also
        :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.

        Parameters :
            bandpass: *str*
                Bandpass for the coefficient you provide.

            coeff: *float*
                Value of the coefficient.

            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts integers
                up to the number of sources, or *None* (i.e., all sources).
                Default is *None*
        """
        if bandpass not in self._bandpasses:
            self._bandpasses.append(bandpass)

        if source is not None:
            self._limb_darkening_coeffs[source - 1].set_limb_coeff_gamma(bandpass, coeff)
        else:
            for i in range(self.n_sources):
                self._limb_darkening_coeffs[i].set_limb_coeff_gamma(bandpass, coeff)

    def get_limb_coeff_gamma(self, bandpass, source=None):
        """
        Get gamma limb darkening coefficient for given band.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient will be provided.

            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts integers
                up to the number of sources, or *None* (i.e., all sources).
                Default is *None*

        Returns :
            gamma: *float* or *list*
                limb darkening coefficient or list of coefficients, in case
                of multiple sources.

        """
        if source is not None and not (1 <= source <= self.n_sources):
            raise ValueError(f'Source number must be between 1 and n_sources = {self.n_sources}.')

        coefficients = self._limb_darkening_coeffs
        if source is not None:
            return coefficients[source - 1].get_limb_coeff_gamma(bandpass)
        elif self.n_sources == 1:
            return coefficients[0].get_limb_coeff_gamma(bandpass)

        return self._get_limb_coeff_gamma_all(bandpass)

    def _get_limb_coeff_gamma_all(self, bandpass):
        """
        Extract limb coeffs for all sources and given band.
        """
        out = []
        for coeff in self._limb_darkening_coeffs:
            try:
                out.append(coeff.get_limb_coeff_gamma(bandpass))
            except KeyError:
                out.append(None)

        if set(out) == set([None]):
            raise ValueError('None of the sources has limb darkening coeff for bandpass ' + bandpass)

        return out

    def _get_limb_coeff_gamma(self, bandpass, gamma, source=None):
        """
        Get gamma from either bandpass or gamma
        """
        if (bandpass is not None) and (gamma is not None):
            raise ValueError('Only one of bandpass and gamma can be set.')
        elif (bandpass is None) and (gamma is None):
            gamma = 0. if self.n_sources == 1 else [0.]*self.n_sources
        elif bandpass is not None:
            if bandpass not in self._bandpasses:
                raise KeyError(f'No limb-darkening coefficient set for {bandpass}.')
            else:
                gamma = self.get_limb_coeff_gamma(bandpass, source)
        else:
            pass

        return gamma

    def set_limb_coeff_u(self, bandpass, coeff, source=None):
        """
        Store u limb darkening coefficient for given band.  See also
        :py:class:`MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient you provide.

            coeff: *float*
                Value of the coefficient.

            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts integers
                up to the number of sources, or *None* (i.e., all sources).
                Default is *None*
        """
        if bandpass not in self._bandpasses:
            self._bandpasses.append(bandpass)

        if source is not None:
            self._limb_darkening_coeffs[source - 1].set_limb_coeff_u(bandpass, coeff)
        else:
            for i in range(self.n_sources):
                self._limb_darkening_coeffs[i].set_limb_coeff_u(bandpass, coeff)

    def get_limb_coeff_u(self, bandpass, source=None):
        """
        Get u limb darkening coefficient for given band.

        Parameters :
            bandpass: *str*
                Bandpass for which coefficient will be provided.

            source: *int* or *None*, optional
                Which source do the given methods apply to? Accepts integers
                up to the number of sources, or *None* (i.e., all sources).
                Default is *None*

        Returns :
            u: *float*
                limb darkening coefficient

        """
        if source is not None and not (1 <= source <= self.n_sources):
            raise ValueError(f'Source number must be between 1 and n_sources = {self.n_sources}.')

        coefficients = self._limb_darkening_coeffs
        if source is not None:
            return coefficients[source - 1].get_limb_coeff_u(bandpass)
        elif self.n_sources == 1:
            return coefficients[0].get_limb_coeff_u(bandpass)
        return [coeff.get_limb_coeff_u(bandpass) for coeff in coefficients]

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
                          bandpass=None, source_flux_ratio=None, separate=None):
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

            source_flux_ratio: *float*, *list*
                If the model has 2 sources, source_flux_ratio is the ratio of
                source_flux_2 / source_flux_1. For N sources, source_flux_ratio a list
                of the ratios of source_flux_i / source_flux_1 where i = 2, N.

            separate: *boolean*, optional
                For multi source models, return magnification of each source
                separately. Defaults to *True* if source_flux_ratio is provided and *False* otherwise
                (then only effective magnification is returned).

        Returns :
            magnification: *np.ndarray*
                A vector of calculated magnification values. For binary source
                models, the effective magnification is returned (unless
                *separate=True*).
        """
        if source_flux_ratio is not None:
            if not isinstance(source_flux_ratio, float):
                raise TypeError(
                    'source_flux_ratio should be a float. Got: {:}'.format(
                        source_flux_ratio))

        gamma = self._get_limb_coeff_gamma(bandpass, gamma)

        if self.n_sources > 1:
            if (source_flux_ratio is None) and (separate is False):
                raise ValueError(
                    'For N sources either source_flux_ratio should be set or' +
                    ' separate=True. \n' +
                    'separate: {0}\n'.format(separate) +
                    'source_flux_ratio: {0}'.format(source_flux_ratio))

            if separate is None:
                if source_flux_ratio is None:
                    separate = True
                else:
                    separate = False

        magnification = self._get_magnification(
            time, satellite_skycoord, gamma, source_flux_ratio, separate)

        return magnification

    def _get_magnification(self, time, satellite_skycoord, gamma,
                           source_flux_ratio, separate):
        """
        Internal function that calculates magnification.
        """
        time = np.atleast_1d(time)

        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(time)

        if gamma is not None and gamma != 0.:
            methods = self.get_magnification_methods()
            self._check_limb_darkening(methods)

        if self.n_sources == 1:
            if source_flux_ratio is not None:
                raise ValueError(
                    'Model.get_magnification() parameter ' +
                    'flux_ratio_constraint has to be None for single source ' +
                    'models, not {:}'.format(source_flux_ratio))
            elif separate is not None:
                raise ValueError(
                    'Model.get_magnification() parameter separate ' +
                    'cannot be True for single source models')
            else:
                magnification = self._magnification_1_source(
                    time, satellite_skycoord, gamma)

        elif self.n_sources >= 2:
            magnification = self._magnification_N_sources(
                time, satellite_skycoord, gamma, source_flux_ratio,
                separate)
        else:
            raise ValueError('Invalid number of sources: {:}'.format(self.n_sources))

        if np.sum(np.isnan(magnification)) > 0:
            fmt = ("EPOCHS:\n{:}\nMODEL:\n{:}Something went wrong with calculating of magnification for " +
                   "the above model. For all above epochs magnifications are NaN.")
            msg = fmt.format(time[np.isnan(magnification)], self.__repr__())
            raise ValueError(msg)

        return magnification

    def get_magnification_curve(self, time, satellite_skycoord, gamma):
        """
        Create a :py:class:`~MulensModel.magnificationcurve.MagnificationCurve`
        object for a given set of times.

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

        Return:
            py:class:`~MulensModel.magnificationcurve.MagnificationCurve`

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

        return magnification_curve

    def _magnification_1_source(self, time, satellite_skycoord, gamma):
        """
        calculate model magnification for given times for model with
        a single source
        """
        magnification_curve = self.get_magnification_curve(
            time, satellite_skycoord, gamma)

        return magnification_curve.get_magnification()

    def _magnification_N_sources(
            self, time, satellite_skycoord, gamma, source_flux_ratio,
            separate):
        """
        calculate model magnification for given times for model with
        two sources

        source_flux_ratio: *float* or *list*
        separate: *bool*
        """
        if separate and (source_flux_ratio is not None):
            raise ValueError(
                'You cannot set both source_flux_ratio and separate' +
                " parameters in Model.get_magnification(). This doesn't " +
                'make sense')

        mags = self._separate_magnifications(time, satellite_skycoord, gamma)

        if separate:
            return mags
        else:
            # Defining source_flux_ratios as relative to source_1 (rather than total flux).
            if isinstance(source_flux_ratio, (float)):
                source_flux_ratio = [source_flux_ratio]

            magnification = mags[0]
            for (mag, flux_ratio) in zip(mags[1:], source_flux_ratio):
                magnification += mag * flux_ratio

            magnification /= (1. + np.sum(source_flux_ratio))
            return magnification

    def get_magnification_curves(self, time, satellite_skycoord, gamma):
        """
        Create a *list* of
        :py:class:`~MulensModel.magnificationcurve.MagnificationCurve`
        objects for multiple sources, given a set of times.

        Parameters :
            time: *np.ndarray*, *list of floats*, or *float*
                Times for which magnification values are requested.

            satellite_skycoord: *astropy.coordinates.SkyCoord*, optional
                *SkyCoord* object that gives satellite positions. Must be
                the same length as time parameter. Use only for satellite
                parallax calculations.

            gamma: *float* or *list*, optional
                The limb-darkening coefficient in gamma convention. Default is 0 which means no limb darkening effect.
                The *list* input allows providing separate gamma for each source.

        Return:
            *list* of
            py:class:`~MulensModel.magnificationcurve.MagnificationCurve`

        """
        kwargs = {'times': time, 'parallax': self._parallax,
                  'coords': self._coords,
                  'satellite_skycoord': satellite_skycoord, 'gamma': gamma}

        mag_curves = []

        for i in range(self.n_sources):
            if isinstance(self._methods, dict):
                methods = self._methods.get(i + 1, None)
            else:
                methods = self._methods

            if isinstance(gamma, (list)):
                kwargs['gamma'] = gamma[i]

            mag_curve = MagnificationCurve(
                 parameters=self.parameters.__getattr__('source_{0}_parameters'.format(i+1)), **kwargs)
            mag_curve.set_magnification_methods(methods, self._default_magnification_method)
            mag_curve.set_magnification_methods_parameters(self._methods_parameters)
            mag_curves.append(mag_curve)

        return mag_curves

    def _separate_magnifications(self, time, satellite_skycoord, gamma):
        """
        Calculate magnification separately for each source.
        """
        mags = []
        mag_curves = self.get_magnification_curves(time, satellite_skycoord, gamma)
        for i in range(self.n_sources):
            self.__setattr__('_magnification_curve_{0}'.format(i + 1), mag_curves[i])
            mag = self.__getattr__('_magnification_curve_{0}'.format(i + 1)).get_magnification()
            mags.append(mag)

        return mags

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
