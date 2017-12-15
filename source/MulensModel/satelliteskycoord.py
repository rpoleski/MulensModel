import numpy as np
from astropy.coordinates import SkyCoord

from MulensModel.horizons import Horizons


class SatelliteSkyCoord(object):
    """
    An object that gives the *Astropy.SkyCoord* of satellite for a given
    epoch based on an ephemerides file.

    Keywords :
        ephemerides_file: *str*
            path to file with satellite ephemerides from JPL horizons,
            for examples see *data/Spitzer_ephemeris_01.dat* or
            *data/K2_ephemeris_01.dat*

        satellite: *str*
            Just the name of the satellite.

    Attributes :
        satellite: *str*
            name of the satellite

    """

    def __init__(self, ephemerides_file, satellite=None):
        """
        ephemerides_file = file with ephemerides for the satellite (Required)
        satellite = Name of the satellite (Optional)
        """
        self.ephemerides_file = ephemerides_file

        self.satellite = satellite

        self._satellite_skycoord = None
        self._horizons = None

    def get_satellite_coords(self, times):
        """
        Calculate the coordinates of the satellite for given times
        using interpolation.

        Parameters :
            times: *np.ndarray* or *list of floats*
                Epochs for which satellite coordinates will be calculated.

        Returns :
            satellite_skycoord: *Astropy.coordinates.SkyCord*
                *SkyCord* for satellite at epochs *times*.

        """
        if self._horizons is None:
            self._horizons = Horizons(self.ephemerides_file)
        # We should add some check if all values in times as within range
        # covered by self._horizons.time.

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
