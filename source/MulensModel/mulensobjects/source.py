from astropy import units as u

from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs


class Source(object):
    """
    Physical properties of a source (background) star.

    Attributes :
        limb_darkening:
        :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`

            Limb darkening coefficients of the source.

    """

    def __init__(self, distance=None, pi_S=None, angular_radius=None,
                 limb_darkening=None):

        self.distance = distance
        self.pi_S = pi_S
        self.angular_radius = angular_radius

        if limb_darkening is None:
            self.limb_darkening = LimbDarkeningCoeffs()
        else:
            self.limb_darkening = limb_darkening

    def __repr__(self):
        return('Source Distance: {0}'.format(self.distance))

    @property
    def distance(self):
        """
        *astropy.Quantity*

        The distance to the source. May be set as a *float*.
        The distance should either be given in pc, or if no unit is
        given, the value is assumed to be kpc (*u.kpc*).

        """
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        if new_distance is None:
            self._distance = new_distance
        else:
            if not isinstance(new_distance, u.Quantity):
                self._distance = new_distance * 1000. * u.pc
            else:
                if (new_distance.unit == "pc") or (new_distance.unit == "kpc"):
                    self._distance = new_distance
                else:
                    raise u.UnitsError(
                        'Allowed units for Source distance are "pc" or "kpc"')

    @property
    def pi_S(self):
        """
        *astropy.Quantity*

        The parallax to the source in millarcseconds. May be set as a
        *float*. If no units are specified, assumes milliarcseconds (*u.mas*).

        """
        return self._distance.to(u.mas, equivalencies=u.parallax())

    @pi_S.setter
    def pi_S(self, new_value):
        if new_value is None:
            pass
        else:
            if not isinstance(new_value, u.Quantity):
                new_value = new_value * u.mas
            self._distance = new_value.to(u.pc, equivalencies=u.parallax())

    @property
    def angular_radius(self):
        """
        *astropy.Quantity*

        Angular radius of the source. May be set as a *float*. If
        units are not specified, assumed to be microarcseconds (*u.mas*).
        """
        return self._angular_radius

    @angular_radius.setter
    def angular_radius(self, new_value):
        if not isinstance(new_value, u.Quantity) and new_value is not None:
            new_value = new_value * u.uas
        self._angular_radius = new_value
