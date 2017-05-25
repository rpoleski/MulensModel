import numpy as np

from astropy.coordinates import SkyCoord#, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.horizons import Horizons

#data_list and ephemrides_file must have the same time standard.
#To implement: mjd2hjd = T/F
#usecols
class MulensData(object):
    """
    A set of photometric measurements for a microlensing event.
    """

    def __init__(self, data_list=None, file_name=None,
                 phot_fmt="mag", coords=None, ra=None, dec=None, 
                 satellite=None, ephemrides_file=None, add_2450000=False,
                 add_2460000=False, bandpass=None, **kwargs):
        """
        Create a MulensData object from a set of photometric measurements.

        Major properties:
           self.time - the dates of the observations
           self.mag - the brightness in magnitudes
           self.err_mag - the errors on the magnitudes
           self.flux - the brightness in flux
           self.err_flux - the errors on the fluxes

        Keywords:
           data_list - a list or array with columns: Date, Magnitude/Flux, Err
           file_name - The path to a file with columns: Date,
               Magnitude/Flux, Err
           *Either data_list or file_name is required.*

           phot_fmt - accepts either 'mag' or 'flux'. Default =
              'mag'. Specifies whether the photometry is in Magnitudes
              or Flux units.

           coords - [optional] sky coordinates of the event
           ra, dec - [optional] sky coordinates of the event
           satellite - [optional] if applicable, specify which
               satellite this dataset comes from.
           ephemrides_file - [optional] specify the ephemrides of the
               satellite when the data were taken. Necessary for
               modeling the satellite parallax effect.

           add_2450000 - Adds 2450000. to the input dates. Useful if
               the dates are supplied as HJD-2450000.
           add_2460000 - Adds 2460000. to the input dates.

           **kwargs - if file_name is provided, uses np.loadtxt to
               load file, and therefore this function accpets loadtxt
               keywords.
        """
        #Initialize some variables
        self._n_epochs = None  
        self._horizons = None
        self._satellite_skycoord = None
        self._init_keys = {'add245':add_2450000, 'add246':add_2460000}
        self._limb_darkening_weights = None
        self.bandpass = bandpass

        #Set the coords (if applicable)...
        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
        #...using coords keyword
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._coords = coords
            else:
                self._coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        #...using ra, dec keywords
        if ra is not None:
            if dec is not None:
                self._coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            else:
                raise AttributeError(coords_msg)
        else:
            if ra is not None:
                raise AttributeError(coords_msg)

        #Import the photometry...
        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            #...from an array
            (vector_1, vector_2, vector_3) = list(data_list) 
            self._initialize(phot_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
        elif file_name is not None:
            #...from a file
            (vector_1, vector_2, vector_3) = np.loadtxt(
                fname=file_name, unpack=True, usecols=(0,1,2), **kwargs)
            self._initialize(phot_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
        
        #Set up satellite properties (if applicable)
        if satellite is None:
            if ephemrides_file is not None:
                raise ValueError(
                    "For datasets with satellite ephemerides file you have"
                    +" to provide satellite name")
            self.is_satellite = False 
        else:
            if ephemrides_file is None:
                raise ValueError(
                    "Currently ephemerides_file has to be specified for each"
                    +" satellite dataset")
            self.ephemrides_file = ephemrides_file
            self.is_satellite = True

    def _initialize(self, phot_fmt, time=None, brightness=None, 
                    err_brightness=None, coords=None):
        """
        internal function to package photometric data using a few numpy arrays

        Keywords:
            phot_fmt - Specifies type of photometry. Either 'flux' or 'mag'. 
            time - Date vector of the data
            brightness - vector of the photometric measurements
            err_brightness - vector of the errors in the phot measurements.
            coords - Sky coordinates of the event
        """
        #Check to see if either add_2450000 or add_2460000 is set 
        #JCY - why is this here? It will create problems once mjd2hjd is
        #implemented, i.e. you might want to do both add_2450000 and mjd2hjd.
        n_additions = 0
        for (key, value) in self._init_keys.items():
            n_additions += value
        if n_additions > 1:
            msg = 'MulensData._initialize(): more than one'
            raise ValueError(msg)

        #Adjust the time vector as necessary.
        if self._init_keys['add245']:
            time += 2450000.
        elif self._init_keys['add246']:
            time += 2460000.

        #Store the time vector
        self._time = time
        self._n_epochs = len(time)

        #Check that the number of epochs equals the number of observations
        if ((len(brightness) != self._n_epochs) 
            or (len(err_brightness) != self._n_epochs)):
            raise ValueError('input data in MulesData have different lengths')

        #Store the photometry
        self._brightness_input = brightness
        self._brightness_input_err = err_brightness        
        self.input_fmt = phot_fmt

        #Create the complementary photometry (mag --> flux, flux --> mag)
        if phot_fmt == "mag":
            self.mag = self._brightness_input
            self.err_mag = self._brightness_input_err
            (self.flux, self.err_flux) = Utils.get_flux_and_err_from_mag(
                                          mag=self.mag, err_mag=self.err_mag)
        elif phot_fmt == "flux":
            self.flux = self._brightness_input
            self.err_flux = self._brightness_input_err
            (self.mag, self.err_mag) = Utils.get_mag_and_err_from_flux(
                                        flux=self.flux, err_flux=self.err_flux)
        else:
            msg = 'unknown brightness format in MulensData'
            raise ValueError(msg)

        #Create an array to flag bad epochs
        self.bad = self.n_epochs * [False]

    @property
    def n_epochs(self):
        """give number of epochs"""
        return self._n_epochs

    @property
    def time(self):
        """short version of time vector"""
        return self._time

    @property
    def coords(self):
        """
        Sky coordinates (RA,Dec)
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        if isinstance(new_value, SkyCoord):
            self._coords = new_value
        else:
            self._coords = SkyCoord(new_value, unit=(u.hourangle, u.deg))

    @property
    def ra(self):
        """
        Right Ascension
        """
        return self._coords.ra

    @ra.setter
    def ra(self, new_value):
        try:
            self._coords.ra = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    new_value, 0.0, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    new_value, self._coords.dec, unit=(u.hourangle, u.deg)) 

    @property
    def dec(self):
        """
        Declination
        """
        return self._coords.dec

    @dec.setter
    def dec(self, new_value):
        try:
            self._coords.dec = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    0.0, new_value, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    self._coords.ra, new_value, unit=(u.hourangle, u.deg))

    @property
    def satellite_skycoord(self):
        """return Astropy SkyCoord of satellite for epochs covered by the dataset"""
        if self.is_satellite is not True:
            #Why not make this return None?
            raise ValueError("You're trying to get satellite information for dataset that has no satellite information")
        if self._satellite_skycoord is None:
            if self._horizons is None:
                self._horizons = Horizons(self.ephemrides_file)
            x = np.interp(
                self._time, self._horizons.time, self._horizons.xyz.x)
            y = np.interp(
                self._time, self._horizons.time, self._horizons.xyz.y)
            z = np.interp(
                self._time, self._horizons.time, self._horizons.xyz.z)
            self._satellite_skycoord = SkyCoord(x=x, y=y, z=z, representation='cartesian')
            self._satellite_skycoord.representation = 'spherical'
        return self._satellite_skycoord

    @property
    def bandpass(self):
        """bandpass of given dataset (primary usage is limb darkening)"""
        return self._bandpass
        
    @bandpass.setter
    def bandpass(self, value):
        if self._limb_darkening_weights is not None:
            raise ValueError("Limb darkening weights were already set - you" +
                                "cannot bandpass now.")
        self._bandpass = value

    def set_limb_darkening_weights(self, weights):
        """save a dictionary weights that will be used to evaluate l
        imb darkening coefficient
        
        e.g. weights = {'I': 1.5, 'V': 1.} if I-band gamma limb-darkening 
        coefficient is 1.5-times larger than V-band"""
        if self.bandpass is not None:
            raise ValueError("Don't try to run " + 
                                "MulensData.set_limb_darkening_weights() " + 
                                "after bandpass was provided")
        if not isinstance(weights, dict):
            raise TypeError("MulensData.set_limb_darkening_weights() " + 
                    "parameter has to be dict, not {:}".format(type(weights)))
        
        self._limb_darkening_weights = weights
        