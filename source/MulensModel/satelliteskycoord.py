import numpy as np
from scipy.interpolate import interp1d
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
        using cubic interpolation. The epheris provided is extrapolated
        up to 3 minutes beyond the range provided.

        Parameters :
            times: *np.ndarray* or *list of floats*
                Epochs for which satellite coordinates will be calculated.

        Returns :
            satellite_skycoord: *Astropy.coordinates.SkyCoord*
                *SkyCoord* for satellite at epochs *times*.

        """
        if self._horizons is None:
            self._prepare_horizons()

        self._check_times(times)

        x = self._interp_x(times)
        y = self._interp_y(times)
        z = self._interp_z(times)

        if int(astropy_version[0]) >= 4:
            self._satellite_skycoord = SkyCoord(
                x=x, y=y, z=z, representation_type='cartesian')
            self._satellite_skycoord.representation_type = 'spherical'
        else:
            self._satellite_skycoord = SkyCoord(
                x=x, y=y, z=z, representation='cartesian')
            self._satellite_skycoord.representation = 'spherical'

        return self._satellite_skycoord

    def _prepare_horizons(self):
        """
        Prepare an instance of Horizons class and interpolation functions.
        """
        self._horizons = Horizons(self._ephemerides_file)
        time = self._horizons.time
        kwargs = dict(kind='cubic', fill_value="extrapolate")
        self._interp_x = interp1d(time, self._horizons.xyz.x, **kwargs)
        self._interp_y = interp1d(time, self._horizons.xyz.y, **kwargs)
        self._interp_z = interp1d(time, self._horizons.xyz.z, **kwargs)

    def _check_times(self, times):
        """
        Make sure that requested range is not too much beyond the ephemeris.
        """
        dt = 3. / (60. * 24)
        min_ = np.min(self._horizons.time)
        max_ = np.max(self._horizons.time)

        if max_ + dt < np.max(times) or min_ - dt > np.min(times):
            msg = ("Satellite ephemeris doesn't cover requested epochs.\n"
                   "Ephemerides file: {:} {:}\nRequested dates: {:} {:}")
            args = [min_, max_, np.min(times), np.max(times)]
            raise ValueError(msg.format(*args))
