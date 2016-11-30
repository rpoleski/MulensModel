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
        pass

    @property
    def distance(self):
        """
        The distance to the lens. An astropy Quantity.
        """
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        """Have not checked what happens if the distance is entered in lyr 
        or AU. Probably we should use new_distance.decompose()"""
        if not isinstance(new_distance, u.Quantity):
            msg1 = 'distance must have astropy units, '
            msg2 = 'i.e. it must be an astropy.units Quantity.'
            raise TypeError(msg1 + msg2)
        else:
            if (new_distance.unit == "pc" or new_distance.unit == "AU" 
                or new_distance.unit == "lyr"):
                self._distance = new_distance
            else:
                raise u.UnitsError('Allowed units are "pc", "AU", or "lyr"') 

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
            
