from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


class Coordinates(SkyCoord):
    """
    A class for the coordinates (RA, Dec) of an event. Inherits from
    astropy.SkyCoord_.

    May be set as a *str*, pair of *str*, or *SkyCoord* object, e.g.

    .. code-block:: python

      Coordinates('18:00:00 -30:00:00')
      Coordinates('18h00m00s', '-30d00m00s')
      Coordinates(SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg)))
      Coordinates(SkyCoord(270.000, -30.000, unit=u.deg))

    If the unit keyword is not specified, defaults to
    unit=(u.hourangle, u.deg) where u is defined by "import
    astropy.units as u".

    .. _astropy.SkyCoord:
      http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html

    """

    def __init__(self,  *args, **kwargs):
        if not isinstance(args[0], SkyCoord) and 'unit' not in kwargs:
            kwargs['unit'] = (u.hourangle, u.deg)
        SkyCoord.__init__(self, *args, **kwargs)
        self._calculate_projected()

    @property
    def galactic_l(self):
        """
        *astropy.coordinates.angles.Longitude*

        Galactic longitude. Note that for connivance, the values l >
        180 degrees are represented as 360-l.
        """
        l = self.galactic.l
        if l > 180. * u.deg:
            l = l - 360. * u.deg
        return l

    @property
    def galactic_b(self):
        """
        *astropy.coordinates.angles.Latitude*

        Galactic latitude calculated from (RA, Dec)
        """
        return self.galactic.b

    @property
    def ecliptic_lon(self):
        """
        *astropy.coordinates.angles.Longitude*

        ecliptic longitude calculated from (RA, Dec)
        """
        from astropy.coordinates import GeocentricTrueEcliptic
        return self.transform_to(GeocentricTrueEcliptic).lon

    @property
    def ecliptic_lat(self):
        """
        *astropy.coordinates.angles.Latitude*

        ecliptic latitude calculated from (RA, Dec)
        """
        from astropy.coordinates import GeocentricTrueEcliptic
        return self.transform_to(GeocentricTrueEcliptic).lat
            
    @property
    def north_projected(self):
        """
        *np.array*

        North projected on the plane of the sky.
        """
        return self._north_projected

    @property
    def east_projected(self):
        """
        *np.array*

        East projected on the plane of the sky.
        """
        return self._east_projected

    def _calculate_projected(self):
        """
        Calculate North and East directions projected on the plane of the sky.
        """
        direction = np.array(self.cartesian.xyz.value)
        north = np.array([0., 0., 1.])
        east_projected = np.cross(north, direction)
        self._east_projected = east_projected / np.linalg.norm(east_projected)
        self._north_projected = np.cross(direction, self._east_projected)

