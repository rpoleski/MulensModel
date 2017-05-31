import numpy as np
from astropy.coordinates import SkyCoord

from MulensModel.horizons import Horizons

class SatelliteSkyCoord(object):
    """
    An object that gives the Astropy SkyCoord of satellite for a given
    epoch based on an ephemerides file.
    """

    def __init__(self, ephemerides_file=None, satellite=None):
        """
        ephemerides_file = file with ephemerides for the satellite (Required)
        satellite = Name of the satellite (Optional)
        """
        if ephemerides_file is None:
            raise ValueError('ephemerides_file must be defined')
        else:
            self.ephemerides_file = ephemerides_file

        self.satellite = satellite
        self._satellite_skycoord = None
        self._horizons = None

    def get_satellite_coords(self, times):
        """
        Return the sky coordinates for the given times by interpolating
        self.ephemerides_file
        """
        if self._horizons is None:
            self._horizons = Horizons(self.ephemerides_file)

        x = np.interp(
                 times, self._horizons.time, self._horizons.xyz.x)
        y = np.interp(
                 times, self._horizons.time, self._horizons.xyz.y)
        z = np.interp(
                 times, self._horizons.time, self._horizons.xyz.z)

        self._satellite_skycoord = SkyCoord(
                  x=x, y=y, z=z, representation='cartesian')
        self._satellite_skycoord.representation = 'spherical'

        return self._satellite_skycoord
