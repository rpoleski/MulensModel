import numpy as np

from astropy import units as u
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

from MulensModel import utils
from MulensModel.modelparameters import ModelParameters
from MulensModel.coordinates import Coordinates
from MulensModel.orbits import Orbit


class Trajectory(object):
    """
    The (dimensionless) X, Y trajectory of the source in the source
    plane. This class includes internal functions that calculate how
    microlensing parallax affects the trajectory.

    For binary lens, the origin of the coordinate system is at the
    center of mass with higher mass (assuming q < 1) at negative X and
    Y=0.

    This class follows the conventions defined in Appendix A of
    `Skowron et al. (2011)
    <https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract>`_
    except the definition of *alpha*, which is shifted by 180 deg.

    Arguments :
        times: [*float*, *list*, *np.ndarray*]
            the times at which to generate the source trajectory,
            e.g. a vector. Either *times* OR (*x* AND *y*) are required.

        parameters: instance of
        :py:class:`~MulensModel.modelparameters.ModelParameters`, required

            a ModelParameters object specifying the microlensing parameters

        x, y: [*float*, *list*, *np.ndarray*]
            Dimensionless X, Y coordinates of trajectory. Either *times* OR
            (*x* AND *y*) are required.

        parallax: *boolean dictionary*, optional
            specifies what parallax effects should be used. Default is
            *False* for each of *'earth_orbital'*, *'satellite'*, and
            *'topocentric'*. (differs from
            :py:class:`~MulensModel.model.Model` which defaults to
            *True*)


        coords: *str*, or
        :py:class:`~MulensModel.coordinates.Coordinates`,
        *Astropy.coordinates.SkyCoord*, optional

            sky coordinates of the event; required for parallax calculations

        satellite_skycoord: *Astropy.coordinates.SkyCoord*, optional
            sky coordinates of the satellite specified by the
            ephemerides file. See
            :py:obj:`~MulensModel.mulensdata.MulensData.satellite_skycoord`.

    Attributes :
        parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
            input :py:class:`~MulensModel.modelparameters.ModelParameters`

        satellite_skycoord: *Astropy.coordinates.SkyCoord*
            input
            :py:attr:`~MulensModel.mulensdata.MulensData.satellite_skycoord`

        parallax: *dict*
            specifies which types of microlensing parallax will be taken
            into account; boolean dict with keys: ``earth_orbital``,
            ``satellite``, and ``topocentric`` (default values are *False*)

        coords: :py:class:`~MulensModel.coordinates.Coordinates`
            event coordinates
    """
    _get_delta_annual_results = dict()
    _get_delta_annual_last = None
    _get_delta_annual_last_index = None
    _get_delta_satellite_results = dict()

    def __init__(self,
                 times=None, parameters=None, x=None, y=None, parallax=None,
                 coords=None, satellite_skycoord=None, earth_coords=None):
        if isinstance(parameters, ModelParameters):
            self.parameters = parameters
        else:
            m = 'parameters is a required and must be a ModelParameters object'
            raise TypeError(m)

        self._set_parallax_and_coords(parallax, coords, satellite_skycoord, earth_coords)

        if times is None:
            if (x is None) and (y is None):
                raise KeyError('Either times or (x AND y) must be defined.')
            self._times = None
            self._x = x
            self._y = y
        else:
            self._times = np.atleast_1d(times)
            self._get_xy()

    def _set_parallax_and_coords(self, parallax, coords, satellite_skycoord, earth_coords):
        """
        Set parallax and coordinates parameters
        """
        self.parallax = {'earth_orbital': False, 'satellite': False, 'topocentric': False}
        if parallax is not None:
            for (key, value) in parallax.items():
                self.parallax[key] = value

        if coords is None or isinstance(coords, Coordinates):
            self.coords = coords
        else:
            self.coords = Coordinates(coords)
        self.satellite_skycoord = satellite_skycoord
        if earth_coords is not None:
            raise NotImplementedError("The earth_coords needed for topocentric parallax is not implemented yet")
        self._earth_coords = None

    @property
    def x(self):
        """
        *np.ndarray*

        Dimensionless X coordinates of trajectory.
        """
        return self._x

    @property
    def y(self):
        """
        *np.ndarray*

        Dimensionless Y coordinates of trajectory.
        """
        return self._y

    @property
    def times(self):
        """
        *np.ndarray* or *None*

        Epochs for which trajectory is calculated.
        *None* if the class was initialized without epochs provided.
        """
        return self._times

    @property
    def parallax_delta_N_E(self):
        """
        *dict*

        Net North (key='N') and East (key='E') components of the parallax
        offset calculated for each time stamp (so sum of the offsets from all
        parallax types).
        """
        return self._delta_N_E

    @property
    def d_perp(self):
        """
        *np.ndarray*

        Projection of the Earth-Satellite
        separation vector (in AU) on the sky.
        """
        if (self.parallax['satellite'] and
                self.satellite_skycoord is not None):
            delta_satellite = self._get_delta_satellite()
            return np.sqrt(delta_satellite['N']**2 + delta_satellite['E']**2)
        else:
            raise AttributeError(
                'Only valid for satellite parallax. ' +
                'satellite_skycoord must be provided.')

    def _get_xy(self):
        """
        For a given set of parameters
        (a :py:class:`~MulensModel.modelparameters.ModelParameters` object),
        calculate the xy position of the source.

        This function has no input and no output. It sets :py:attr:`~x` and
        :py:attr:`~y` attributes.
        """
        # Calculate the position of the source
        vector_tau = (self._times - self.parameters.t_0) / self.parameters.t_E
        vector_u = self.parameters.u_0 * np.ones(self._times.size)

        if 'pi_E_N' in self.parameters.parameters:
            shifts = self._get_shifts_parallax()
            vector_tau += shifts[0]
            vector_u += shifts[1]

        if self.parameters.is_xallarap:
            shifts = self._get_shifts_xallarap()
            vector_tau += shifts[0]
            vector_u += shifts[1]

        # If 2 lenses, rotate trajectory relative to binary lens axis
        is_mass_sheet = self.parameters.is_external_mass_sheet_with_shear
        n_lenses = self.parameters.n_lenses
        if n_lenses == 1 and not is_mass_sheet:
            vector_x = vector_tau
            vector_y = vector_u
        elif n_lenses == 2 or (n_lenses == 1 and is_mass_sheet):
            alpha = self.parameters.get_alpha(self._times) * (np.pi / 180)
            sin_alpha = np.sin(alpha)
            cos_alpha = np.cos(alpha)

            vector_x = -vector_tau * cos_alpha + vector_u * sin_alpha
            vector_y = -vector_tau * sin_alpha - vector_u * cos_alpha
            # The above equations use alpha in counterclockwise
            # convention, i.e., the same as proposed by Skowron et
            # al. (2011), but shifted by 180 deg.
        else:
            raise NotImplementedError(
                "trajectory for more than 2 lenses not handled yet")

        self._x = vector_x
        self._y = vector_y

    def _get_shifts_parallax(self):
        """calculate shifts caused by parallax effect"""
        if self.coords is None:
            raise ValueError(
                "You're trying to calculate trajectory in a parallax model, "
                "but event sky coordinates were not provided.")

        keys = ['earth_orbital', 'satellite', 'topocentric']
        if set([self.parallax[k] for k in keys]) == set([False]):
            raise ValueError(
                'If pi_E value is provided then at least one value of '
                'parallax dict has to be True '
                '(earth_orbital, satellite, or topocentric)')

        self._calculate_delta_N_E()
        return self._project_delta()

    def _calculate_delta_N_E(self):
        """
        Calculate shifts caused by microlensing parallax effect.
        """
        self._delta_N_E = {'N': 0., 'E': 0.}

        if self.parallax['earth_orbital']:
            delta_annual = self._get_delta_annual()
            self._delta_N_E['N'] += delta_annual['N']
            self._delta_N_E['E'] += delta_annual['E']

        if (self.parallax['satellite'] and
                self.satellite_skycoord is not None):
            delta_satellite = self._get_delta_satellite()
            self._delta_N_E['N'] += delta_satellite['N']
            self._delta_N_E['E'] += delta_satellite['E']

        if self.parallax['topocentric'] and self._earth_coords is not None:
            # When you implement it, make sure the behavior depends on the
            # access to the observatory location information as the satellite
            # parallax depends on the access to satellite_skycoord.
            raise NotImplementedError(
                "The topocentric parallax effect not implemented yet")

    def _project_delta(self):
        """
        Project N and E parallax offset vector onto the tau, beta plane.
        """
        delta = self.parallax_delta_N_E

        delta_tau = (delta['N'] * self.parameters.pi_E_N +
                     delta['E'] * self.parameters.pi_E_E)
        delta_beta = (-delta['N'] * self.parameters.pi_E_E +
                      delta['E'] * self.parameters.pi_E_N)
        return [delta_tau, delta_beta]

    def _get_delta_annual(self):
        """
        calculates projected Earth positions required by annual parallax
        """
        index = (self.parameters.t_0_par, self.coords.ra.value,
                 self.coords.dec.value, tuple(self._times.tolist()))
        if index == Trajectory._get_delta_annual_last_index:
            return Trajectory._get_delta_annual_last
        if index in Trajectory._get_delta_annual_results:
            return Trajectory._get_delta_annual_results[index]
        time_ref = self.parameters.t_0_par

        velocity = utils.Utils.velocity_of_Earth(time_ref) / 1731.45683
        # We change units from km/s to AU/d.

        if not np.all(np.isfinite(self._times)):
            msg = "Some times have incorrect values: {:}".format(self._times[~np.isfinite(self._times)])
            raise ValueError(msg)

        position = get_body_barycentric(
            body='earth', time=Time(self._times, format='jd', scale='tdb'))
        position_ref = get_body_barycentric(
            body='earth', time=Time(time_ref, format='jd', scale='tdb'))
        # Seems that get_body_barycentric depends on time system, but there is
        # no way to set BJD part of BJD_TDB in astropy.Time(). The option
        # *format* above indicates if the first argument of Time() is
        # a float indicating JD or e.g., a string in the form
        # '1999-01-01T00:00:00.123' - this would be value 'fits'.
        # Hence, the user has to provide BJD times (or at least HJD).

        # Main calculation is in 2 lines below:
        delta_s = (position_ref.xyz.T - position.xyz.T).to(u.au).value
        delta_s += np.outer(self._times - time_ref, velocity)
        # and the results require projecting on the plane of the sky:
        out_n = np.dot(delta_s, self.coords.north_projected)
        out_e = np.dot(delta_s, self.coords.east_projected)

        out = {'N': out_n, 'E': out_e}
        Trajectory._get_delta_annual_results[index] = out
        Trajectory._get_delta_annual_last_index = index + tuple()
        Trajectory._get_delta_annual_last = out
        return out

    def _get_delta_satellite(self):
        """
        calculates differences of Earth and satellite positions
        projected on the plane of the sky at event position
        """
        index = (self.coords.ra.value, self.coords.dec.value,
                 tuple(self._times.tolist()))
        if index in Trajectory._get_delta_satellite_results.keys():
            return Trajectory._get_delta_satellite_results[index]

        # Transform the satellite ephemerides.
        satellite = self.satellite_skycoord
        satellite.transform_to(frame=self.coords.frame)

        # Project the satellite parallax effect based on the direction of
        # the event.
        direction = np.array(self.coords.cartesian.xyz.value)
        north_projected = self.coords.north_projected
        east_projected = self.coords.east_projected
        delta_satellite = {}
        dot = utils.Utils.dot
        delta_satellite['N'] = -dot(satellite.cartesian, north_projected).value
        delta_satellite['E'] = -dot(satellite.cartesian, east_projected).value
        delta_satellite['D'] = -dot(satellite.cartesian, direction).value

        Trajectory._get_delta_satellite_results[index] = delta_satellite
        return delta_satellite

    def _get_shifts_xallarap(self):
        """calculate shifts caused by xallarap effect"""
        zip_ = self.parameters.parameters.items()
        t_0_xi = self.parameters.t_0_xi
        orbit_parameters = {
            key[3:]: value for (key, value) in zip_ if key[:3] == "xi_"}
        orbit_parameters['epoch_reference'] = t_0_xi
        orbit = Orbit(**orbit_parameters)
        positions = orbit.get_reference_plane_position(self._times)
        return positions - self.parameters.xallarap_reference_position
