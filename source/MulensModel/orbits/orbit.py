import numpy as np
import math

"""
Classes that can be accessed directly:
- Orbit
- OrbitCircular
- OrbitEccentric
"""


class Orbit(object):
    """
    Class that combines :py:class:`OrbitCircular` and
    :py:class:`OrbitEccentric`.
    """
    def __new__(self, **kwargs):
        for Class_ in [OrbitCircular, OrbitEccentric]:
            try:
                new = Class_(**kwargs)
            except Exception:
                pass
            else:
                return new
        raise RuntimeError("Orbit.__new__() failed")


class _OrbitAbstract(object):
    """
    Abstract class for orbits.
    """
    def _check_circular_orbit_parameters(self, semimajor_axis):
        """
        Check if period and semimajor axis make physical sense.
        """
        if self._period <= 0.:
            raise ValueError('Orbital period has to be positive.\n'
                             'Provided value: ' + str(self._period))
        if semimajor_axis <= 0.:
            raise ValueError('Semimajor axis has to be positive.\n'
                             'Provided value: ' + str(semimajor_axis))

    def _check_for_and_get_periapsis_epoch(
            self, periapsis_epoch, argument_of_latitude_reference,
            epoch_reference):
        """
        Check if arguments properly define epoch of
        the periapsis (eccentric orbit) or
        the ascending node (circular orbit) passage
        """
        if periapsis_epoch is not None:
            if argument_of_latitude_reference is not None:
                raise RuntimeError(
                    "periapsis_epoch and argument_of_latitude_reference "
                    "cannot be both set")
            if epoch_reference is not None:
                raise RuntimeError(
                    "periapsis_epoch and epoch_reference cannot be both set")
            return periapsis_epoch

        if argument_of_latitude_reference is None or epoch_reference is None:
            raise RuntimeError("Not enough arguments to define the epoch of "
                               "periapsis/ascending node")

        u_reference = argument_of_latitude_reference * np.pi / 180.
        return self._get_periapsis_epoch(u_reference, epoch_reference)

    def _set_circular_orbit_parameters(self, period, semimajor_axis,
                                       Omega_node, inclination,
                                       periapsis_epoch):
        """
        Set parameters that are used for circular orbits.
        """
        self._semimajor_axis = semimajor_axis
        Omega = math.pi * Omega_node / 180.
        self._rotation_matrix = np.array([[math.cos(Omega), -math.sin(Omega)],
                                          [math.sin(Omega), math.cos(Omega)]])
        self._cos_inclination = math.cos(math.pi * inclination / 180.)
        self._periapsis_epoch = periapsis_epoch

    def get_reference_plane_position(self, time):
        """
        Calculate positions for given epochs in the reference plane.

        Parameters :
            time: *float* or *np.ndarray*
                Epochs for which positions are requested.

        Returns :
            positions: *np.ndarray*
                Positions of the body at given epochs.
        """
        projected_position = self._get_projected_position(time)
        position = np.matmul(self._rotation_matrix, projected_position)
        return position

    def _get_projected_position(self, time):
        """
        Get position projected on reference plane,
        but not rotated to the reference coordinate system
        """
        position = self.get_orbital_plane_position(time)
        position[1, ] *= self._cos_inclination
        return position

    def _get_eccentric_anomaly(self, time):
        """
        Calculate eccentric anomaly (typically indicated by E).
        """
        mean_anomaly = self._get_mean_anomaly(time)
        anomaly = self._get_normalized_anomaly_minus_pi_pi(mean_anomaly)
        eccentric_anomaly = (
            self._get_eccentric_anomaly_from_normalized_mean_anomaly(anomaly))
        return eccentric_anomaly

    def _get_mean_anomaly(self, time):
        """
        Calculate mean anomaly, i.e., the one that is linear in time and
        typically indicated by M.
        """
        anomaly = 2. * math.pi * (time - self._periapsis_epoch) / self._period
        return anomaly

    def _get_normalized_anomaly_minus_pi_pi(self, anomaly):
        """
        get value normalized to (-pi, pi) range
        """
        out = self._get_normalized_anomaly_zero_two_pi(anomaly + np.pi) - np.pi
        return out

    def _get_normalized_anomaly_zero_two_pi(self, anomaly):
        """
        get value normalized to (0, 2*pi) range
        """
        return np.remainder(anomaly, 2 * np.pi)

    def _get_eccentric_anomaly_from_normalized_mean_anomaly(self,
                                                            mean_anomaly):
        """
        Turn mean anomaly in range (-pi, pi) into eccentric anomaly
        """
        anomaly = mean_anomaly + 0.
        if self._eccentricity > 0.95:  # This limit was manually found.
            anomaly = np.pi * np.sign(mean_anomaly)
        for _ in range(5):
            anomaly += self._get_anomaly_correction(anomaly, mean_anomaly)
        return anomaly

    def _get_anomaly_correction(self, anomaly, mean_anomaly):
        """
        Calculations needed to solve Kepler's equation.
        We use Newton's method.
        The input anomaly is current estimate of eccentric anomaly.
        """
        numerator = (
            anomaly - self._eccentricity * np.sin(anomaly) - mean_anomaly)
        denominator = 1. - self._eccentricity * np.cos(anomaly)
        return -numerator / denominator


class OrbitCircular(_OrbitAbstract):
    """
    Circular orbit.

    Keywords :
        period: *float*
            Orbital period of binary in days

        semimajor_axis: *float*
            Semimajor axis of the orbit. The unit is not specified.
            Note that the positions returned by
            :py:func:`get_orbital_plane_position()` and
            :py:func:`get_reference_plane_position()`
            functions will be in the same units.

        Omega_node: *float*
            Longitude of the ascending node, i.e., the angle from
            the reference direction to the ascending node direction.

        inclination: *float*
            Inclination of the orbit relative to plane of the sky.

        ascending_node_epoch: *float* or *None*
            Epoch when body is in the ascending node.
            It's in days and usually you want to provide full BJD or HJD.

        argument_of_latitude_reference: *float* or *None*
            Argument of latitude (i.e., u = omega + nu(t_ref)) for
            *epoch_reference*, which together define
            *ascending_node_epoch* (omega).

        epoch_reference: *float* or *None*
            Reference epoch that together with
            *argument_of_latitude_reference* defines
            *ascending_node_epoch* (omega).
    """
    def __init__(self, period, semimajor_axis, Omega_node, inclination,
                 ascending_node_epoch=None,
                 argument_of_latitude_reference=None, epoch_reference=None):
        self._period = period
        self._check_circular_orbit_parameters(semimajor_axis)
        ascending_node_epoch = self._check_for_and_get_periapsis_epoch(
            ascending_node_epoch, argument_of_latitude_reference,
            epoch_reference)
        self._set_circular_orbit_parameters(
            period, semimajor_axis, Omega_node,
            inclination, ascending_node_epoch)

    def _get_periapsis_epoch(self, u_reference, epoch_reference):
        """
        Calculate ascending node epoch
        (called periapis in general case of eccentric orbits) based on
        the argument_of_latitude (u) at given epoch
        """
        time_shift = self._period * u_reference / (2. * np.pi)
        return np.float64(epoch_reference - time_shift)

    def get_orbital_plane_position(self, time):
        """
        Calculate positions in the orbital plane for given epochs

        Parameters :
            time: *float* or *np.ndarray*
                Epochs for which positions are requested.

        Returns :
            positions: *np.ndarray*
                Calculated positions.
        """
        anomaly = self._get_mean_anomaly(time)
        unit_vector = np.array([np.cos(anomaly), np.sin(anomaly)])
        return self._semimajor_axis * unit_vector


class OrbitEccentric(_OrbitAbstract):
    """
    Eccentric orbit.

    Keywords :
        period: *float*
            Orbital period of binary.

        semimajor_axis: *float*
            Semimajor axis of the orbit. The unit is not specified.
            Note that the positions returned by
            :py:func:`get_orbital_plane_position()` and
            :py:func:`get_reference_plane_position()`
            functions will be in the same units.

        Omega_node: *float*
            Longitude of the ascending node, i.e., the angle from
            the reference direction to the ascending node direction.

        inclination: *float*
            Inclination of the orbit relative to plane of the sky.

        eccentricity: *float*
            Eccentricity of the orbit, has to be in (0, 1) range.

        omega_periapsis: *float*
            Argument of periapsis in degrees.

        periapsis_epoch: *float* or *None*
            Epoch when body is in periapsis.
            It's in days and usually you want to provide full BJD or HJD.

        argument_of_latitude_reference: *float* or *None*
            Argument of latitude (i.e., u = omega + nu(t_ref)) for
            *epoch_reference*, which together define
            *periapsis_epoch* (omega).

        epoch_reference: *float* or *None*
            Reference epoch that together with
            *argument_of_latitude_reference* defines
            *periapsis_epoch* (omega).

    """
    def __init__(
            self, period, semimajor_axis, Omega_node, inclination,
            eccentricity, omega_periapsis, periapsis_epoch=None,
            argument_of_latitude_reference=None, epoch_reference=None):
        self._period = period
        self._omega_periapsis = omega_periapsis * np.pi / 180.
        self._eccentricity = eccentricity
        self._check_circular_orbit_parameters(semimajor_axis)
        periapsis_epoch = self._check_for_and_get_periapsis_epoch(
            periapsis_epoch, argument_of_latitude_reference, epoch_reference)
        self._set_circular_orbit_parameters(
            period, semimajor_axis, Omega_node, inclination, periapsis_epoch)

    def _get_periapsis_epoch(self, u_reference, epoch_reference):
        """
        Calculate periapsis epoch (omega) based on
        the argument_of_latitude (u) at given reference epoch
        """
        true_anomaly = u_reference - self._omega_periapsis
        mean_anomaly = self._get_mean_anomaly_from_true_anomaly(true_anomaly)
        mean_anomaly = self._get_normalized_anomaly_minus_pi_pi(mean_anomaly)
        time_shift = self._period * mean_anomaly / (2. * np.pi)
        return epoch_reference - time_shift

    def _get_mean_anomaly_from_true_anomaly(self, true_anomaly):
        """
        Calculate mean anomaly (M) based on true anomaly (E)
        """
        (sin_E, cos_E) = self._get_sin_cos_eccentric_anomaly(true_anomaly)
        E = np.arctan2(sin_E, cos_E)
        return E - self._eccentricity * sin_E

    def _get_sin_cos_eccentric_anomaly(self, true_anomaly):
        """
        Calculate sin(E) and cos(E) based on true_anomaly (nu)
        """
        cos_nu = np.cos(true_anomaly)
        denominator = 1 + self._eccentricity * cos_nu
        sin_E = (np.sqrt(1 - self._eccentricity**2) * np.sin(true_anomaly)
                 / denominator)
        cos_E = (self._eccentricity + cos_nu) / denominator
        return (sin_E, cos_E)

    def get_orbital_plane_position(self, time):
        """
        Calculate positions in the orbital plane for given epochs

        Parameters :
            time: *float* or *np.ndarray*
                Epochs for which positions are requested.

        Returns :
            positions: *np.ndarray*
                Calculated positions.
        """
        eccentric_anomaly = self._get_eccentric_anomaly(time)
        out_x = np.cos(eccentric_anomaly) - self._eccentricity
        out_y = np.sqrt(1 - self._eccentricity**2) * np.sin(eccentric_anomaly)
        return self._semimajor_axis * np.array([out_x, out_y])

    def get_true_anomaly_deg(self, time):
        """
        Calculate true anomaly [deg] for given epochs.

        Parameteres :
            time: *float* or *np.ndarray*
                Epochs for which positions are requested.

        Returns :
            true_anomaly: *float* or *np.ndarray*
                Values of true anomaly (nu) for given epochs.
                The results are in 0-360 range.
        """
        true_anomaly = self._get_true_anomaly(time)
        true_anomaly_wrapped = np.remainder(true_anomaly, 2*np.pi)
        return true_anomaly_wrapped * (180. / np.pi)

    def _get_true_anomaly(self, time):
        """
        Calculate true anomaly for given times.
        The result is in (-pi, pi) range.
        """
        eccentric_anomaly = self._get_eccentric_anomaly(time)
        (sin_nu, cos_nu) = self._get_sin_cos_true_anomaly(eccentric_anomaly)
        return np.arctan2(sin_nu, cos_nu)

    def _get_sin_cos_true_anomaly(self, eccentric_anomaly):
        """
        Calculate sin(nu) and cos(nu) based on E
        """
        cos_E = np.cos(eccentric_anomaly)
        denominator = 1. - self._eccentricity * cos_E
        cos_nu = (cos_E - self._eccentricity) / denominator
        sin_nu = (np.sqrt(1-self._eccentricity**2) * np.sin(eccentric_anomaly)
                  / denominator)
        return (sin_nu, cos_nu)
