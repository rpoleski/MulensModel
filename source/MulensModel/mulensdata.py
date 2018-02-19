import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.coordinates import Coordinates


# data_list and ephemerides_file must have the same time standard.
# To implement: mjd2hjd = T/F
# usecols
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

    Keywords :
        data_list: [*list* of *lists*, *numpy.ndarray*], optional
            columns: Date, Magnitude/Flux, Err

        file_name: *str*, optional
            The path to a file with columns: Date, Magnitude/Flux,
            Err. Loaded using np.loadtxt. See ``**kwargs``.

        **Either data_list or file_name is required.**

        phot_fmt: *str*
           Specifies whether the photometry is in Magnitudes or Flux
           units. accepts either 'mag' or 'flux'. Default = 'mag'.

        coords: *astropy.SkyCoord*, optional
           sky coordinates of the event

        ra, dec: *str*, optional
           sky coordinates of the event

        ephemerides_file: *str*, optional
           Specify the ephemerides of a satellite over the period when
           the data were taken. Will be interpolated as necessary to
           model the satellite parallax effect. See "Instructions on
           getting satellite positions" in MulensModel.README.md

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
            Flags for good data, should be the same length as the number of 
            data points.

        ``**kwargs`` - :py:func:`np.loadtxt()` keywords. Used if
        file_name is provided.

    Attributes (all vectors):

        flux - the brightness in flux

        err_flux - the errors on the fluxes

    """

    def __init__(self, data_list=None, file_name=None,
                 phot_fmt="mag", coords=None, ra=None, dec=None,
                 ephemerides_file=None, add_2450000=False,
                 add_2460000=False, bandpass=None, bad=None, good=None,
                 **kwargs):

        # Initialize some variables
        self._n_epochs = None
        self._horizons = None
        self._satellite_skycoord = None

        self._init_keys = {'add245': add_2450000, 'add246': add_2460000}
        self._limb_darkening_weights = None
        self.bandpass = bandpass

        # Set the coords (if applicable)...
        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
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

        # Import the photometry...
        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            # ...from an array
            (vector_1, vector_2, vector_3) = list(data_list)
            self._initialize(
                phot_fmt, time=vector_1, brightness=vector_2,
                err_brightness=vector_3, coords=self._coords)
        elif file_name is not None:
            # ...from a file
            (vector_1, vector_2, vector_3) = np.loadtxt(
                fname=file_name, unpack=True, usecols=(0, 1, 2), **kwargs)
            self._initialize(
                phot_fmt, time=vector_1, brightness=vector_2,
                err_brightness=vector_3, coords=self._coords)
        else:
            raise ValueError(
                'MulensData cannot be initialized with ' +
                'data_list or file_name')

        if bad is not None and good is not None:
            raise ValueError('Provide bad or good, but not both')
        elif bad is not None:
            self.bad = bad
        elif good is not None:
            self.good = good
        else:
            self.bad = self.n_epochs * [False]

        # Set up satellite properties (if applicable)
        self.ephemerides_file = ephemerides_file

    def _initialize(self, phot_fmt, time=None, brightness=None,
                    err_brightness=None, coords=None):
        """
        Internal function to import photometric data into the correct
        form using a few numpy ndarrays.

        Parameters:
            phot_fmt - Specifies type of photometry. Either 'flux' or 'mag'.
            time - Date vector of the data
            brightness - vector of the photometric measurements
            err_brightness - vector of the errors in the phot measurements.
            coords - Sky coordinates of the event, optional

        """
        if self._init_keys['add245'] and self._init_keys['add246']:
            raise ValueError(
                'You cannot initialize MulensData with both ' +
                'add_2450000 and add_2460000 being True')

        # Adjust the time vector as necessary.
        if self._init_keys['add245']:
            time += 2450000.
        elif self._init_keys['add246']:
            time += 2460000.

        # Store the time vector
        if time.dtype != np.float64:
            raise TypeError((
                    'time vector in MulensData() must be of ' +
                    'numpy.float64 type, not {:}').format(time.dtype))
        self._time = time
        self._n_epochs = len(time)

        # Check that the number of epochs equals the number of observations
        if ((len(brightness) != self._n_epochs) or
                (len(err_brightness) != self._n_epochs)):
            raise ValueError('input data in MulesData have different lengths')

        # Store the photometry
        self._brightness_input = brightness
        self._brightness_input_err = err_brightness
        self.input_fmt = phot_fmt

        # Create the complementary photometry (mag --> flux, flux --> mag)
        if phot_fmt == "mag":
            self._mag = self._brightness_input
            self._err_mag = self._brightness_input_err
            (self.flux, self.err_flux) = Utils.get_flux_and_err_from_mag(
                mag=self.mag, err_mag=self.err_mag)
        elif phot_fmt == "flux":
            self.flux = self._brightness_input
            self.err_flux = self._brightness_input_err
            self._mag = None
            self._err_mag = None
        else:
            msg = 'unknown brightness format in MulensData'
            raise ValueError(msg)

    @property
    def bad(self):
        """
        *np.ndarray boolean*

        flags marking bad data
        """
        return self._bad
        
    @bad.setter
    def bad(self, new_value):
        self._bad = new_value
        self._good = np.logical_not(self._bad)
        
    @property
    def good(self):
        """
        *np.ndarray boolean*

        flags marking good data i.e., opposite to py:func:`bad`
        """
        return self._good

    @good.setter
    def good(self, new_value):
        self._good = new_value
        self._bad = np.logical_not(self._good)

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
    def coords(self):
        """
        see :py:class:`~MulensModel.coordinates.Coordinates`
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        self._coords = Coordinates(new_value)

    @property
    def n_epochs(self):
        """
        *int*

        give total number of epochs (including bad data)
        """
        return self._n_epochs

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
    def bandpass(self):
        """
        *string*

        bandpass of given dataset (primary usage is limb darkening), e.g. 'I'
        or 'V'
        """
        return self._bandpass

    @bandpass.setter
    def bandpass(self, value):
        if self._limb_darkening_weights is not None:
            raise ValueError(
                "Limb darkening weights were already set - you" +
                "cannot bandpass now.")
        self._bandpass = value

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

