import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, exists
import warnings

from MulensModel.utils import Utils, PlotUtils
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.coordinates import Coordinates


class MulensData(object):
    """
    A set of photometric measurements for a microlensing event.

    Examples of how to define a MulensData object:
        data = MulensData(file_name=SAMPLE_FILE_01)

        data = MulensData(data_list=[[Dates], [Magnitudes], [Errors]])

    **Parallax calculations assume that the dates supplied are
    BJD_TDB. See** :py:class:`~MulensModel.trajectory.Trajectory`. If
    you aren't using parallax, the time system shouldn't matter as
    long as it is consistent across all MulensData and Model objects.
    If you have multiple datasets, then you also need multiple instances
    of MulensData class.

    Keywords :
        data_list: [*list* of *lists*, *numpy.ndarray*], optional
            The list that contains three *lists* or *numpy.ndarrays*
            that specify: time, magnitude or flux, and its uncertainty
            (in that order). The lengths of these three objects must be
            the same.

        file_name: *str*, optional
            The path to a file with columns: Date, Magnitude/Flux,
            Err. Loaded using :py:func:`numpy.loadtxt()`. See ``**kwargs``.

        **Either data_list or file_name is required.**

        phot_fmt: *str*
           Specifies whether the photometry is provided in magnitude or flux
           space. Accepts either 'mag' or 'flux'. Default = 'mag'.

        chi2_fmt: *str*
           Specifies whether the format used for chi^2 calculation
           should be done in Magnitude or Flux spaces. Accepts either
           'mag' or 'flux'. Default is 'flux' because almost always
           the errors are gaussian in flux space.

        coords: *astropy.SkyCoord*, optional
           sky coordinates of the event

        ra, dec: *str*, optional
           sky coordinates of the event

        ephemerides_file: *str*, optional
           Specify the ephemerides of a satellite over the period when
           the data were taken. You may want to extend the time range
           to get nicer plots. Will be interpolated as necessary to
           model the satellite parallax effect. See instructions_ on
           getting satellite positions.
           Note that there is no check on time format (e.g., BJD TBD vs. HJD)
           and it should be the same as in *data_list* or *file_name*.

        add_2450000: *boolean*, optional
            Adds 2450000 to the input dates. Useful if the dates
            are supplied as HJD-2450000.

        add_2460000: *boolean*, optional
            Adds 2460000 to the input dates. Useful if the dates
            are supplied as HJD-2460000.

        bandpass: see :obj:`bandpass`

        bad: *boolean np.ndarray*, optional
            Flags for bad data (data to exclude from fitting and
            plotting). Should be the same length as the number of data
            points.

        good: *boolean np.ndarray*, optional
            Flags for good data, should be the same length as the
            number of data points.

        plot_properties: *dict*, optional
            Specify properties for plotting, e.g. ``color``, ``marker``,
            ``label``, ``alpha``, ``zorder``, ``markersize``, ``visible``,
            and also the ``show_bad`` and ``show_errorbars``
            properties.

            Note: pyplot functions errorbar() and scatter() are used to
            plot data with errorbars and without them, respectively.
            The type and size of marker are specified using different
            keywords: ('fmt', 'markersize') for errorbar() and
            ('marker', 'size') for scatter(). You can use either convention
            in :py:attr:`plot_properties` and they will be translated
            to appropriate keywords. If there are similar problems with
            other keywords, then they won't be translated unless you
            contact code authors.

            Other special keys :
                show_errorbars: *boolean*, optional
                    Whether or not to show the errorbars for this dataset.

                show_bad: *boolean*, optional
                    Whether or not to plot data points flagged as bad.

        ``**kwargs``:
            Kwargs passed to np.loadtxt(). Works only if ``file_name`` is set.

    .. _instructions:
        https://github.com/rpoleski/MulensModel/blob/master/documents/Horizons_manual.md

    """

    def __init__(self, data_list=None, file_name=None,
                 phot_fmt="mag", chi2_fmt="flux",
                 coords=None, ra=None, dec=None,
                 ephemerides_file=None, add_2450000=False,
                 add_2460000=False, bandpass=None, bad=None, good=None,
                 plot_properties=None, **kwargs):

        self._n_epochs = None
        self._horizons = None
        self._satellite_skycoord = None

        self._init_keys = {'add245': add_2450000, 'add246': add_2460000}
        self._limb_darkening_weights = None
        self.bandpass = bandpass
        self._chi2_fmt = chi2_fmt
        self._file_name = file_name
        self._input_fmt = phot_fmt

        self._set_coords(coords=coords, ra=ra, dec=dec)

        if plot_properties is None:
            plot_properties = dict()
        self._plot_properties = plot_properties

        self._import_photometry(data_list, **kwargs)

        if bad is not None and good is not None:
            raise ValueError('Provide bad or good, but not both')
        elif bad is not None:
            self.bad = bad
        elif good is not None:
            self.good = good
        else:
            self.bad = self.n_epochs * [False]

        # Set up satellite properties (if applicable)
        self._ephemerides_file = ephemerides_file

    def _import_photometry(self, data_list, **kwargs):
        """import time, brightnes, and its uncertainy"""
        # Import the photometry...
        if data_list is not None and self._file_name is not None:
            raise ValueError(
                'MulensData cannot be initialized with both data_list and ' +
                'file_name. Choose one or the other.')
        elif data_list is not None:  # ...from an array
            if len(kwargs) > 0:
                raise ValueError('data_list and kwargs cannot be both set')
            if len(data_list) != 3:
                try:
                    msg0 = "\n" + str(data_list) + "\n"
                except Exception:
                    msg0 = ""
                msg = (msg0 + "\n" +
                       'MulensData was initiated with data_list of length ' +
                       '{:}, while length of 3 is expected (i.e. time, mag ' +
                       'or flux, and uncertainty).')
                raise ValueError(msg.format(len(data_list)))
            (vector_1, vector_2, vector_3) = list(data_list)
            self._initialize(
                time=np.array(vector_1), brightness=np.array(vector_2),
                err_brightness=np.array(vector_3), coords=self._coords)
        elif self._file_name is not None:  # ...from a file
            usecols = kwargs.pop('usecols', (0, 1, 2))
            if not exists(self._file_name):
                raise FileNotFoundError(self._file_name)
            try:
                (vector_1, vector_2, vector_3) = np.loadtxt(
                    fname=self._file_name, unpack=True,
                    usecols=usecols, **kwargs)
            except Exception:
                print("kwargs passed to np.loadtxt():")
                print(kwargs)
                print("usecols =", usecols)
                print("File:", self._file_name)
                raise
            self._initialize(
                time=vector_1, brightness=vector_2,
                err_brightness=vector_3, coords=self._coords)

            # Check if data label specified, if not use file_name.
            if 'label' not in self.plot_properties.keys():
                if self._file_name is not None:
                    self.plot_properties['label'] = basename(self._file_name)
                else:
                    self.plot_properties['label'] = 'a dataset'
        else:
            raise ValueError(
                'MulensData cannot be initialized with ' +
                'data_list or file_name')

    def _initialize(self, time=None, brightness=None,
                    err_brightness=None, coords=None):
        """
        Internal function to import photometric data into the correct
        form using a few numpy ndarrays.

        Parameters:
            time - Date vector of the data
            brightness - vector of the photometric measurements
            err_brightness - vector of the errors in the phot measurements.
            coords - Sky coordinates of the event, optional

        """
        if self._init_keys['add245'] and self._init_keys['add246']:
            raise ValueError(
                'You cannot initialize MulensData with both ' +
                'add_2450000 and add_2460000 being True')

        if time.dtype != np.float64:
            raise TypeError((
                'time vector in MulensData() must be of ' +
                'numpy.float64 type, not {:}').format(time.dtype))

        # Adjust the time vector as necessary.
        if self._init_keys['add245']:
            time += 2450000.
        elif self._init_keys['add246']:
            time += 2460000.

        # Store the time vector
        self._time = time
        self._n_epochs = len(time)

        # Check that the number of epochs equals the number of observations
        if ((len(brightness) != self._n_epochs) or
                (len(err_brightness) != self._n_epochs)):
            raise ValueError('input data in MulesData have different lengths')

        # Store the photometry
        self._brightness_input = brightness
        self._brightness_input_err = err_brightness

        # Create the complementary photometry (mag --> flux, flux --> mag)
        if self._input_fmt == "mag":
            self._mag = self._brightness_input
            self._err_mag = self._brightness_input_err
            (self._flux, self._err_flux) = Utils.get_flux_and_err_from_mag(
                mag=self.mag, err_mag=self.err_mag)
            if np.min(self._err_flux) <= 0.:
                msg = ("Scaling of magnitude uncertainties to flux space "
                       "resulted in zero or negative values. Maybe the "
                       "photometry format is in fact 'flux', not 'mag' "
                       "(as you indicated). ")
                if self._file_name is not None:
                    msg += "File name: " + self._file_name
                raise ValueError(msg)
        elif self._input_fmt == "flux":
            self._flux = self._brightness_input
            self._err_flux = self._brightness_input_err
            self._mag = None
            self._err_mag = None
        else:
            msg = 'unknown brightness format in MulensData'
            raise ValueError(msg)

    def _set_coords(self, coords=None, ra=None, dec=None):
        """Set the coordinates and raise errors if applicable."""
        self._coords = None
        if (coords is not None) or (ra is not None) or (dec is not None):
            # Check for errors and if none, set the coordinates
            warnings.warn(
                'coords will be deprecated in future. There is no reason ' +
                'to tie this to a given dataset', FutureWarning)
            coords_msg = 'Must specify both or neither of ra and dec'
            # ...using coords keyword
            if coords is not None:
                self._coords = Coordinates(coords)
            # ...using ra, dec keywords
            if ra is not None:
                if dec is not None:
                    self._coords = Coordinates(ra, dec)
                else:
                    raise AttributeError(coords_msg)
            else:
                if ra is not None:
                    raise AttributeError(coords_msg)

    @property
    def plot_properties(self):
        """
        *dict*

        Settings that specify how the photometry should be plotted.

        The keys in this *dict* could be either special keys introduced here
        (i.e., ``show_bad`` and ``show_errorbars``) or keys accepted by
        matplotlib.pyplot plotting functions. The latter could be for example
        ``color``, ``marker``, ``label``, ``alpha``, ``zorder``,
        ``markersize``, or ``visible``.

        See :py:class:`~MulensModel.mulensdata.MulensData`
        for more information.
        """
        return self._plot_properties

    def plot(self, phot_fmt=None, show_errorbars=None, show_bad=None,
             subtract_2450000=False, subtract_2460000=False,
             model=None, plot_residuals=False, **kwargs):
        """
        Plot the data.

        Uses :py:attr:`plot_properties` for label, color, etc.
        This settings can be changed by setting ``**kwargs``.

        You can plot in either flux or magnitude space.

        Keywords:
            phot_fmt: *string* ('mag', 'flux')
                Whether to plot the data in magnitudes or in flux. Default
                is the same as :py:attr:`input_fmt`.

            show_errorbars: *boolean*
                If show_errorbars is True (default), plots with
                matplotlib.errorbar(). If False, plots with
                matplotlib.scatter().

            show_bad: *boolean*
                If False, bad data are suppressed (default).
                If True, shows points marked as bad
                (:py:obj:`mulensdata.MulensData.bad`) as 'x'

            subtract_2450000, subtract_2460000: *boolean*
                If True, subtracts 2450000 or 2460000 from the time
                axis to get more human-scale numbers. If using it, make
                sure to also set the same settings for all other
                plotting calls (e.g. :py:func:`plot_lc()`).

            model: :py:class:`~MulensModel.model.Model`
                DEPRECATED. Use :py:func:`~MulensModel.model.Event.plot_data()`
                to plot a dataset scaled to a model.

            plot_residuals: *boolean*
                If *True* then residuals are plotted (*model* is required).
                Default is *False*, i.e., plot the data.

            ``**kwargs``:
                passed to matplotlib plotting functions.
        """
        if phot_fmt is None:
            phot_fmt = self.input_fmt
        if phot_fmt not in ['mag', 'flux']:
            raise ValueError('wrong value of phot_fmt: {:}'.format(phot_fmt))
        if plot_residuals and model is None:
            raise ValueError(
                'MulensData.plot() requires model to plot residuals')

        if model is None:
            (y_value, y_err) = self._get_y_value_y_err(phot_fmt,
                                                       self.flux,
                                                       self.err_flux)
        else:
            raise KeyError(
                'Passing a model to MulensData.plot will be depracated. Use ' +
                'Event.plot_data() or Event.plot_residuals() instead.')

        self._plot_datapoints(
            (y_value, y_err), subtract_2450000=subtract_2450000,
            subtract_2460000=subtract_2460000, show_errorbars=show_errorbars,
            show_bad=show_bad, **kwargs)

        if phot_fmt == 'mag':
            (ymin, ymax) = plt.gca().get_ylim()
            if ymax > ymin:
                plt.gca().invert_yaxis()

    def _plot_datapoints(
            self, y, subtract_2450000=False,
            subtract_2460000=False, show_errorbars=None, show_bad=None,
            **kwargs):
        """
        plot datapoints while evaluating various contingencies
        """
        (y_value, y_err) = y
        subtract = PlotUtils.find_subtract(subtract_2450000=subtract_2450000,
                                           subtract_2460000=subtract_2460000)

        if show_errorbars is None:
            show_errorbars = self.plot_properties.get('show_errorbars', True)

        if show_bad is None:
            show_bad = self.plot_properties.get('show_bad', False)

        properties = self._set_plot_properties(
            show_errorbars=show_errorbars, **kwargs)
        properties_bad = self._set_plot_properties(
            show_errorbars=show_errorbars, bad=True, **kwargs)
        if 'label' in properties_bad.keys():
            properties_bad['label'] = None

        time_good = self.time[self.good] - subtract
        time_bad = self.time[self.bad] - subtract

        if show_errorbars:
            container = self._plt_errorbar(time_good, y_value[self.good],
                                           y_err[self.good], properties)
            if show_bad:
                if 'color' in properties_bad or 'c' in properties_bad:
                    pass
                else:
                    properties_bad['color'] = container[0].get_color()

                self._plt_errorbar(time_bad, y_value[self.bad],
                                   y_err[self.bad], properties_bad)
        else:
            collection = self._plt_scatter(time_good, y_value[self.good],
                                           properties)
            if show_bad:
                change = True
                keys = ['c', 'color', 'facecolor', 'facecolors', 'edgecolors']
                for key in keys:
                    change &= key not in properties_bad
                if change:
                    properties_bad['color'] = collection.get_edgecolor()
                self._plt_scatter(time_bad, y_value[self.bad], properties_bad)

    def _set_plot_properties(self, show_errorbars=True, bad=False, **kwargs):
        """
        Set plot properties using ``**kwargs`` and
        `py:plot_properties`. kwargs takes precedent.

        Keywords:
            show_errorbars: *boolean*
                `True` means plotting done with plt.errorbar. `False`
                means plotting done with plt.scatter.

            bad: *boolean*
                `True` means marker is default to 'x'. `False` means
                marker is default to 'o'.

           ``**kwargs``: *dict*
               Keywords accepted by plt.errorbar() or plt.scatter().
        """
        if show_errorbars:
            marker_key = 'fmt'
            size_key = 'markersize'  # In plt.errorbar(), 'ms' is equivalent.
        else:
            marker_key = 'marker'
            size_key = 's'
        marker_keys_all = ['marker', 'fmt']
        size_keys_all = ['markersize', 'ms', 's']

        # Some older versions of matplotlib have problems when both
        # 'fmt' and 'color' are specified. Below we take a list of formats
        # from Notes section of:
        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html
        if 'fmt' in kwargs:
            for char in ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']:
                if char in kwargs['fmt']:
                    kwargs['fmt'] = kwargs['fmt'].replace(char, "")
                    kwargs['color'] = char

        properties = {}

        # Overwrite dataset settings (i.e., self.plot_properties) with kwargs.
        for dictionary in [self.plot_properties, kwargs]:
            for (key, value) in dictionary.items():
                if key in marker_keys_all:
                    properties[marker_key] = value
                elif key in size_keys_all:
                    properties[size_key] = value
                else:
                    properties[key] = value

        if bad:
            properties[marker_key] = 'x'
        elif marker_key not in properties.keys():
            properties[marker_key] = 'o'

        if size_key not in properties.keys():
            properties[size_key] = 5

        for remove_key in ['show_bad', 'show_errorbars']:
            properties.pop(remove_key, None)

        return properties

    def _plt_errorbar(self, time, y, yerr, kwargs):
        """
        save run of matplotlib.pyplot.errorbar(); returns ErrorbarContainer
        """
        try:
            container = plt.errorbar(time, y, yerr=yerr, **kwargs)
        except Exception:
            print("kwargs passed to plt.errorbar():")
            print(kwargs)
            raise
        return container

    def _plt_scatter(self, time, y, kwargs):
        """
        save run of matplotlib.pyplot.scatter(); returns PathCollection
        """
        try:
            collection = plt.scatter(time, y, **kwargs)
        except Exception:
            print("kwargs passed to plt.scatter():")
            print(kwargs)
            raise
        return collection

    def _get_y_value_y_err(self, phot_fmt, flux, flux_err):
        """
        just calculate magnitudes if needed, or return input otherwise
        """
        if phot_fmt == 'mag':
            return Utils.get_mag_and_err_from_flux(flux, flux_err)
        else:
            return (flux, flux_err)

    def set_limb_darkening_weights(self, weights):
        """
        Save a dictionary of weights that will be used to evaluate the
        limb darkening coefficient. See also
        :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`

        Parameters :
            weights: *dict*
                A dictionary that specifies weight for each
                bandpass. Keys are *str* and values are *float*, e.g.,
                ``{'I': 1.5, 'V': 1.}`` if the I-band gamma
                limb-darkening coefficient is 1.5-times larger than
                the V-band.

        """
        if self.bandpass is not None:
            raise ValueError(
                "Don't try to run MMulensData.set_limb_darkening_weights() " +
                "after bandpass was provided")
        if not isinstance(weights, dict):
            raise TypeError(
                "MulensData.set_limb_darkening_weights() " +
                "parameter has to be dict, not {:}".format(type(weights)))

        self._limb_darkening_weights = weights

    @property
    def coords(self):
        """
        :py:class:`~MulensModel.coordinates.Coordinates`

        Sky coordinates of data.
        See :py:class:`~MulensModel.coordinates.Coordinates`.
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        self._coords = Coordinates(new_value)

    @property
    def time(self):
        """
        *np.ndarray*

        vector of dates
        """
        return self._time

    @property
    def mag(self):
        """
        *np.ndarray*

        magnitude vector
        """
        if self._mag is None:
            (self._mag, self._err_mag) = Utils.get_mag_and_err_from_flux(
                flux=self.flux, err_flux=self.err_flux)
        return self._mag

    @property
    def err_mag(self):
        """
        *np.ndarray*

        vector of magnitude errors
        """
        if self._err_mag is None:
            self.mag
        return self._err_mag

    @property
    def flux(self):
        """
        *numpy.ndarray*

        Vector of the measured brightness in flux units.
        """
        if self._flux is None:
            (self._flux, self._err_flux) = Utils.get_flux_and_err_from_mag(
                mag=self.mag, err_mag=self.err_mag)
        return self._flux

    @property
    def err_flux(self):
        """
        *np.ndarray*

        Vector of uncertainties of *flux* values.
        """
        if self._err_flux is None:
            self.flux

        return self._err_flux

    @property
    def bad(self):
        """
        *np.ndarray boolean*

        flags marking bad data
        """
        return self._bad

    @bad.setter
    def bad(self, new_value):
        new_value = np.asarray(new_value)
        if new_value.dtype != np.dtype('bool'):
            raise TypeError("MulensData.bad has to be a boolean numpy array")

        self._bad = new_value
        self._good = np.logical_not(self._bad)

    @property
    def good(self):
        """
        *np.ndarray boolean*

        flags marking good data i.e., opposite to :py:func:`bad`
        """
        return self._good

    @good.setter
    def good(self, new_value):
        new_value = np.asarray(new_value)
        if new_value.dtype != np.dtype('bool'):
            raise TypeError("MulensData.good has to be a boolean numpy array")

        self._good = new_value
        self._bad = np.logical_not(self._good)

    @property
    def n_epochs(self):
        """
        *int*

        give total number of epochs (including bad data)
        """
        return self._n_epochs

    def data_and_err_in_input_fmt(self):
        """
        Gives photometry in input format (mag or flux).

        Returns :
            data: *np.ndarray*
                Magnitudes or fluxes

            data_err: *np.ndarray*
                Uncertainties of magnitudes or of fluxes

        """
        return self._get_data_and_err_in_fmt(self.input_fmt)

    def _get_data_and_err_in_fmt(self, fmt):
        """
        get data and their photometry in mag or flux
        """
        if fmt == "mag":
            data = self.mag
            err_data = self.err_mag
        elif fmt == "flux":
            data = self.flux
            err_data = self.err_flux
        else:
            raise ValueError('Unrecognized data format: {:}'.format(fmt))
        return (data, err_data)

    def data_and_err_in_chi2_fmt(self):
        """
        Gives photometry in format used for chi2 calculation
        (flux in most cases, but magnitude possible).

        Returns :
            data: *np.ndarray*
                Magnitudes or fluxes

            data_err: *np.ndarray*
                Uncertainties of magnitudes or of fluxes

        """
        return self._get_data_and_err_in_fmt(self.chi2_fmt)

    @property
    def bandpass(self):
        """
        *String*

        Bandpass of given dataset (primary usage is limb darkening), e.g. 'I'
        or 'V'. Returns *None* if not set.
        """
        return self._bandpass

    @bandpass.setter
    def bandpass(self, value):
        if self._limb_darkening_weights is not None:
            raise ValueError(
                "Limb darkening weights were already set - you" +
                "cannot bandpass now.")

        self._bandpass = value

    @property
    def satellite_skycoord(self):
        """
        *Astropy.SkyCoord* object for satellite
        positions at epochs covered by the dataset

        Returns :
            skycoord: *astropy.coordinates.SkyCoord*
                satellite positions at epochs covered by the dataset
        """
        if self.ephemerides_file is None:
            raise ValueError('ephemerides_file is not defined.')

        if self._satellite_skycoord is None:
            satellite_skycoord = SatelliteSkyCoord(
                ephemerides_file=self.ephemerides_file)
            self._satellite_skycoord = satellite_skycoord.get_satellite_coords(
                self._time)

        return self._satellite_skycoord

    @property
    def input_fmt(self):
        """
        *str* ('mag' or 'flux')

        Input format - same as *phot_fmt* keyword in __init__().
        """
        return self._input_fmt

    @property
    def chi2_fmt(self):
        """
        *str* ('mag' or 'flux')

        Photometry format used  for chi^2 calculations. Default is 'flux'.
        """
        return self._chi2_fmt

    @property
    def ephemerides_file(self):
        """
        *str*

        File with satellite ephemeris.
        """
        return self._ephemerides_file

    def copy(self):
        """
        Returns a copy of given instance with settings copied

        Returns :
            mulens_data: :py:class:`~MulensModel.mulensdata.MulensData`
                Copy of self.
        """
        data_and_err = self.data_and_err_in_input_fmt()
        kwargs = {
            'data_list': [self.time, *list(data_and_err)],
            'phot_fmt': self.input_fmt, 'chi2_fmt': self._chi2_fmt,
            'coords': self.coords, 'ephemerides_file': self._ephemerides_file,
            'add_2450000': self._init_keys['add245'],
            'add_2460000': self._init_keys['add246'],
            'bandpass': self.bandpass, 'bad': np.array(self.bad),
            'plot_properties': {**self.plot_properties}
            }

        out = MulensData(**kwargs)
        out._file_name = self._file_name

        return out
