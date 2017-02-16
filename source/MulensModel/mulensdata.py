import sys
import numpy as np

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.horizons import Horizons

#data_list and ephemrides_file must have the same time standard.
#To implement: mjd2hjd = T/F
class MulensData(object):
    def __init__(self, data_list=None, file_name=None,
                 mag_fmt="mag", coords=None, ra=None, dec=None, 
                 satellite=None, ephemrides_file=None, add_2450000=False,
                 add_2460000=False):
        self._n_epochs = None  
        self._horizons = None
        self._satellite_skycoord = None
        self._init_keys = {'add245':add_2450000, 'add246':add_2460000}

        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._coords = coords
            else:
                self._coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        if ra is not None:
            if dec is not None:
                self._coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            else:
                raise AttributeError(coords_msg)
        else:
            if ra is not None:
                raise AttributeError(coords_msg)

        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            vector_1, vector_2, vector_3 = list(data_list) 
            self._initialize(mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
        elif file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(
                fname=file_name, unpack=True, usecols=(0,1,2))
            self._initialize(mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
        
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

    def _initialize(self, mag_fmt, time=None, brightness=None, 
                    err_brightness=None, coords=None):
        """internal function to initialized data using a few numpy arrays"""
        n_additions = 0
        for (key, value) in self._init_keys.items():
            n_additions += value
        if n_additions > 1:
            msg = 'More than one delta time found in MulensData._initialize()'
            raise ValueError(msg)

        if self._init_keys['add245']:
            time += 2450000.
        elif self._init_keys['add246']:
            time += 2460000.

        self._time = time
        self._n_epochs = len(time)

        if ((len(brightness) != self._n_epochs) 
            or (len(err_brightness) != self._n_epochs)):
            raise ValueError('input data in MulesData have different lengths')

        self._brightness_input = brightness
        self._brightness_input_err = err_brightness        
        self.input_fmt = mag_fmt

        if mag_fmt == "mag":
            self.mag = self._brightness_input
            self.err_mag = self._brightness_input_err
            (self.flux, self.err_flux) = Utils.get_flux_and_err_from_mag(
                                          mag=self.mag, err_mag=self.err_mag)
        elif mag_fmt == "flux":
            self.flux = self._brightness_input
            self.err_flux = self._brightness_input_err
            (self.mag, self.err_mag) = Utils.get_mag_and_err_from_flux(
                                        flux=self.flux, err_flux=self.err_flux)
        else:
            msg = 'unknown format of brightness in ' + file_name + ' file'
            raise ValueError(msg)

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

