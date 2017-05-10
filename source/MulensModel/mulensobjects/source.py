from astropy import units as u

from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs

class Source(object):
    """
    Physical properties of a source (background) star.
    """
    def __init__(self, distance=None, angular_size=None,
                 limb_darkening=None):
        self.distance=distance
        self.angular_size=angular_size
        if limb_darkening is None:
            self.limb_darkening = LimbDarkeningCoeffs()

    def __repr__(self):
        return('Source Distance: {0}'.format(self.distance))

    @property
    def distance(self):
        """
        The distance to the lens. An astropy Quantity.
        """
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        """
        The distance should either be given in pc, or if no unit is
        given, the value is assumed to be kpc if it is <50 and in pc
        otherwise.
        """
        if new_distance is None:
            self._distance = new_distance
        else:
            if not isinstance(new_distance, u.Quantity):
                if new_distance < 50:
                    self._distance = new_distance * 1000. * u.pc
                else:
                    self._distance = new_distance * u.pc
            else:
                if (new_distance.unit == "pc") or (new_distance.unit == "kpc"):
                    self._distance = new_distance
                else:
                    raise u.UnitsError(
                        'Allowed units for Source distance are "pc" or "kpc"') 

    @property
    def pi_S(self):
        """
        The parallax to the lens in millarcseconds.
        """
        return self._distance.to(u.mas, equivalencies=u.parallax())

    @pi_S.setter
    def pi_S(self, new_value):
        if not isinstance(new_value, u.Quantity):
            new_value = new_value * u.mas
        self._distance = new_value.to(u.pc, equivalencies=u.parallax())
            
        
