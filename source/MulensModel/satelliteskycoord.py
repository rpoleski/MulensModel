import numpy as np
from scipy.interpolate import interp1d
import warnings
from astropy.coordinates import SkyCoord
from astropy import __version__ as astropy_version

from MulensModel.horizons import Horizons


class SatelliteSkyCoord(object):
    """
    An object that gives the *Astropy.SkyCoord* of satellite for a given
    epoch based on an ephemerides file.

    Keywords :
        ephemerides_file: *str*
            path to file with satellite ephemerides from JPL horizons,
            for examples see *data/ephemeris_files/Spitzer_ephemeris_01.dat*
            or *data/ephemeris_files/K2_ephemeris_01.dat*

        satellite: *str*, optional
            Just the name of the satellite.

    Attributes :
        satellite: *str*
            name of the satellite

    """

    def __init__(self, ephemerides_file, satellite=None):
        self._ephemerides_file = ephemerides_file

        self.satellite = satellite

        self._satellite_skycoord = None
        self._horizons = None

    def get_satellite_coords(self, times):
        """
        Calculate the coordinates of the satellite for given times
        using cubic interpolation.

        Parameters :
            times: *np.ndarray* or *list of floats*
                Epochs for which satellite coordinates will be calculated.

        Returns :
            satellite_skycoord: *Astropy.coordinates.SkyCoord*
                *SkyCoord* for satellite at epochs *times*.

        """
        if self._horizons is None:
            self._horizons = Horizons(self._ephemerides_file)

        time = self._horizons.time
        if (np.max(time) + 0.001 < np.max(times) or
                np.min(time) - 0.001 > np.min(times)):
            msg_1 = "Ephemerides file: {:} {:}\n ".format(
                np.min(time), np.max(time))
            msg_2 = "Requested dates: {:} {:}".format(
                np.min(times), np.max(times))
            raise ValueError(
                "Satellite ephemeris doesn't cover requested epochs.\n " +
                msg_1 + msg_2)

        x = interp1d(time, self._horizons.xyz.x, kind='cubic')(times)
        y = interp1d(time, self._horizons.xyz.y, kind='cubic')(times)
        z = interp1d(time, self._horizons.xyz.z, kind='cubic')(times)

        if int(astropy_version[0]) >= 4:
            self._satellite_skycoord = SkyCoord(
                x=x, y=y, z=z, representation_type='cartesian')
            self._satellite_skycoord.representation_type = 'spherical'
        else:
            self._satellite_skycoord = SkyCoord(
                x=x, y=y, z=z, representation='cartesian')
            self._satellite_skycoord.representation = 'spherical'

        return self._satellite_skycoord
