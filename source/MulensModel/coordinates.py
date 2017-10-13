from astropy.coordinates import SkyCoord
from astropy import units as u

class Coordinates(SkyCoord):
    """
    A class for Sky Coordinates.
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

if __name__ == "__main__":
    coords = Coordinates("18:00:00 -30:00:00", unit=(u.hourangle, u.deg))
    print(coords)

    coords_1 = SkyCoord("18:00:00 -30:00:00", unit=(u.hourangle, u.deg))
    coords_2 = Coordinates(coords_1)
    print(coords_2)

    coords_3 = Coordinates("17:00:00 -28:00:00")
    print(coords_3)

    coords_4 = Coordinates("17:00:00", "-27:00:00")
    print(coords_4)
