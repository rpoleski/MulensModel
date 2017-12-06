from astropy.coordinates import SkyCoord
from astropy import units as u

class Coordinates(SkyCoord):
    """
    A class for the coordinates (RA, Dec) of an event. Inherits from
    astropy.SkyCoord_.

    May be set as a *str*, pair of *str*, or *SkyCoord* object, e.g.

        Coordinates('18:00:00 -30:00:00')

        Coordinates('18h00m00s', '-30d00m00s')

        Coordinates(SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg)))
        
    If the unit keyward is not specified, defaults to unit=(u.hourangle, u.deg)

    where u is defined by "import astropy.units as u". 

    .. _astropy.SkyCoord: http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    """

    def __init__(self,  *args, **kwargs):
        if not isinstance(args[0], SkyCoord) and not 'unit' in kwargs:
            kwargs['unit'] = (u.hourangle, u.deg)
        SkyCoord.__init__(self, *args, **kwargs)

    @property
    def galactic_l(self):
        """Galactic longitude. Note that for connivance, the values l > 180 
        degrees are represented as 360-l."""
        l = self.galactic.l
        if l > 180. * u.deg:
            l = l - 360. * u.deg
        return l

    @property
    def galactic_b(self):
        """Galactic latitude calculated from (RA, Dec)"""
        return self.galactic.b

    @property
    def ecliptic_lon(self):
        """ecliptic longitude calculated from (RA, Dec)"""
        from astropy.coordinates import GeocentricTrueEcliptic
        return self.transform_to(GeocentricTrueEcliptic).lon

    @property
    def ecliptic_lat(self):
        """ecliptic latitude calculated from (RA, Dec) """
        from astropy.coordinates import GeocentricTrueEcliptic
        return self.transform_to(GeocentricTrueEcliptic).lat
