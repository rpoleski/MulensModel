import numpy as np
import warnings

from astropy.coordinates import SkyCoord, ICRS, FK4, FK5
from astropy.coordinates import GeocentricTrueEcliptic
from astropy import units as u

from MulensModel.utils import Utils


class Coordinates(SkyCoord):
    """
    A class for the coordinates (RA, Dec) of an event. Inherits from
    astropy.SkyCoord_.

    May be set as a *str*, pair of *str*, or *SkyCoord* object, e.g.

    .. code-block:: python

      from astropy.coordinates import SkyCoord
      from astropy import units as u
      Coordinates('18:00:00 -30:00:00')
      Coordinates('18h00m00s', '-30d00m00s')
      Coordinates(SkyCoord('18:00:00 -30:00:00',
                  unit=(u.hourangle, u.deg)))
      Coordinates(SkyCoord(270.000, -30.000, unit=u.deg))

    If the unit keyword is not specified, defaults to
    unit=(u.hourangle, u.deg) where u is defined by "import
    astropy.units as u".

    You can print an instance of this class.

    .. _astropy.SkyCoord:
      http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html

    """

    def __init__(self,  *args, **kwargs):
        if isinstance(args[0], str):
            if 'unit' not in kwargs and len(args) > 0:
                self._check_for_ra_in_degrees(args[0])
                kwargs['unit'] = (u.hourangle, u.deg)
        self._validate_input(args[0], kwargs.get('frame'))

        SkyCoord.__init__(self, *args, **kwargs)
        if self.cartesian.xyz.shape not in [(3,), (3, 1)]:
            raise ValueError(
                "Something wrong with parameters of Coordinates().\nMost " +
                "probably you have provided more than one sky position.")
        self._calculate_projected()

    def _check_for_ra_in_degrees(self, arg):
        """
        Try to check if RA is outside 0-24 range and raise a warning if that
        is the case. No warning would be given if the user provided RA of,
        e.g., 5.123 and assumed degrees, but has not explicitly stated that.
        """
        if isinstance(arg, str):
            arg = arg.split()[0]

        try:
            value = float(arg)
        except Exception:
            return

        if value > 24.:
            warning = (
                'It seems that you provided RA in degrees rather than hours. '
                'We suggest to provide RA unit or use hours as RA units. ' +
                str(value))
            warnings.warn(warning, UserWarning)
        elif value < 0.:
            warning = (
                "It's very uncommon to use negative RA. Please remember that "
                "a default unit for RA is hours (not degrees). " + str(value))
            warnings.warn(warning, UserWarning)

    def _validate_input(self, arg, frame):
        """
        Validate input for coordinates, checking if format is allowed (ICRS,
        FK4 or FK5) or raising ValueError otherwise. If SkyCoord() instance
        is provided, the frame should be allowed as well.
        """
        allowed_fmts = (str, SkyCoord, ICRS, FK4, FK5)
        if not isinstance(arg, allowed_fmts):
            raise ValueError('Coordinate format {:} is not allowed.'.format(type(arg)))

        if isinstance(arg, str):
            if frame not in [None, 'icrs', 'fk4', 'fk5']:
                raise ValueError("Only ICRS, FK4 and FK5 frames are allowed" +
                                 " to Coordinates().")

        elif isinstance(arg, (SkyCoord, u.quantity.Quantity)):
            test = '18h00m00s -30d00m00s'
            is_icrs = arg.is_equivalent_frame(SkyCoord(test, frame='icrs'))
            is_fk4 = arg.is_equivalent_frame(SkyCoord(test, frame='fk4'))
            is_fk5 = arg.is_equivalent_frame(SkyCoord(test, frame='fk5'))
            if not (is_icrs | is_fk4 | is_fk5):
                raise ValueError("Provided SkyCoord is not in allowed frame.")

    def _calculate_projected(self):
        """
        Calculate North and East directions projected on the plane of the sky.
        """
        direction = np.array(self.cartesian.xyz.value)
        if direction.shape == (3, 1):
            direction = direction[:, 0]
        north = np.array([0., 0., 1.])
        self._east_projected = Utils.vector_product_normalized(
            north, direction)
        self._north_projected = np.cross(direction, self._east_projected)

    def __repr__(self):
        return self.to_string(style="hmsdms", sep=":", precision=2)

    @property
    def galactic_l(self):
        """
        *astropy.coordinates.angles.Longitude*

        Galactic longitude. Note that for connivance, the values l >
        180 degrees are represented as 360-l.
        """
        gal_l = SkyCoord(self).galactic.l
        if gal_l > 180. * u.deg:
            gal_l = gal_l - 360. * u.deg
        return gal_l

    @property
    def galactic_b(self):
        """
        *astropy.coordinates.angles.Latitude*

        Galactic latitude calculated from (RA, Dec)
        """
        return SkyCoord(self).galactic.b

    @property
    def ecliptic_lon(self):
        """
        *astropy.coordinates.angles.Longitude*

        ecliptic longitude calculated from (RA, Dec)
        """
        return SkyCoord(self).transform_to(GeocentricTrueEcliptic).lon

    @property
    def ecliptic_lat(self):
        """
        *astropy.coordinates.angles.Latitude*

        ecliptic latitude calculated from (RA, Dec)
        """
        return SkyCoord(self).transform_to(GeocentricTrueEcliptic).lat

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

    def v_Earth_projected(self, full_BJD):
        """
        Earth velocity at *full_BJD* projected on the plane of sky towards
        given coordinates.

        Parameters :
            full_BJD: *float*
                Epoch for which projected velocity is requested. In most cases
                it is
                :py:attr:`~MulensModel.modelparameters.ModelParameters.t_0_par`

        Returns :
            v_Earth_perp_N: *float*
                North component of Earth's projected velocity in km/s.

            v_Earth_perp_E: *float*
                East component of Earth's projected velocity in km/s.
        """
        velocity = Utils.velocity_of_Earth(full_BJD)

        v_Earth_perp_N = np.dot(velocity, self.north_projected)
        v_Earth_perp_E = np.dot(velocity, self.east_projected)

        return (v_Earth_perp_N, v_Earth_perp_E)
