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
        for Class_ in [OrbitCircular, OrbitEccentric, OrbitEccentricThieleInnes]:
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
        self._Omega_node = Omega
        self._rotation_matrix_reference = np.array([[math.cos(Omega), -math.sin(Omega)],
                                                    [math.sin(Omega), math.cos(Omega)]])
        self._inclination = math.pi * inclination / 180.
        self._cos_inclination = math.cos(self._inclination)
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
        position = np.matmul(self._rotation_matrix_reference, projected_position)
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
        self._rotation_matrix_orbital = np.array([[math.cos(self._omega_periapsis), -math.sin(self._omega_periapsis)],
                                                  [math.sin(self._omega_periapsis), math.cos(self._omega_periapsis)]])
        self._argument_of_latitude_reference = argument_of_latitude_reference
        self._eccentricity = eccentricity
        self._check_circular_orbit_parameters(semimajor_axis)
        self._epoch_reference = epoch_reference
        self._argument_of_latitude_reference = argument_of_latitude_reference
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
        # The ksi is towards pericenter and eta is perpendicular to it.
        ksi = np.cos(eccentric_anomaly) - self._eccentricity
        eta = np.sqrt(1 - self._eccentricity**2) * np.sin(eccentric_anomaly)
        orbital_plane_versor = np.matmul(self._rotation_matrix_orbital, np.array([ksi, eta]))
        return self._semimajor_axis * orbital_plane_versor

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

    def thiele_innes_orbit_elements_dict_degrees(self):
        """
        Return standard orbital elements and Thiele-Innes constants in a dictionary with angles in degrees.
        """
        sin_Omega_node = math.sin(self._Omega_node)
        cos_Omega_node = math.cos(self._Omega_node)
        sin_omega_periapsis = math.sin(self._omega_periapsis)
        cos_omega_periapsis = math.cos(self._omega_periapsis)
        cos_inclination = self._cos_inclination
        A = cos_Omega_node * cos_omega_periapsis - sin_Omega_node * sin_omega_periapsis * cos_inclination
        B = sin_Omega_node * cos_omega_periapsis + cos_Omega_node * sin_omega_periapsis * cos_inclination
        F = -cos_Omega_node * sin_omega_periapsis - sin_Omega_node * cos_omega_periapsis * cos_inclination
        G = -sin_Omega_node * sin_omega_periapsis + cos_Omega_node * cos_omega_periapsis * cos_inclination
        A *= self._semimajor_axis
        B *= self._semimajor_axis
        F *= self._semimajor_axis
        G *= self._semimajor_axis
        dict_out = self.orbit_elements_dict_degrees()
        dict_out['A'] = A
        dict_out['B'] = B
        dict_out['F'] = F
        dict_out['G'] = G
        return dict_out

    def orbit_elements_dict_degrees(self):
        """
        Return standard orbital elements in a dictionary with angles in degrees.
        """
        dict_out = {'period': self._period, 'semimajor_axis': self._semimajor_axis, 'eccentricity': self._eccentricity,
                    'inclination': self._inclination * 180. / np.pi,
                    'omega_periapsis': self._omega_periapsis * 180. / np.pi,
                    'Omega_node': self._Omega_node * 180. / np.pi}
        if self._argument_of_latitude_reference is not None:
            dict_out['argument_of_latitude_reference'] = self._argument_of_latitude_reference
            dict_out['epoch_reference'] = self._epoch_reference
        dict_out['periapsis_epoch'] = self._periapsis_epoch

        return dict_out


class OrbitEccentricThieleInnes(OrbitEccentric):
    """
    Class for eccentric orbits defined by Thiele-Innes constants.
    Based on `An Introduction to Close Binary Stars, R. W. Hilditch (2001)`

    Keywords:
        period: *float*
            Orbital period of binary.

        semimajor_axis: *float*
            Semimajor axis of the orbit. The unit is not specified.
            Note that the positions returned by
            : py: func: `get_orbital_plane_position()` and
            : py: func: `get_reference_plane_position()`
            functions will be in the same units.

        eccentricity: *float*
            Eccentricity of the orbit, has to be in (0, 1) range.

        A: *float*
            Thiele-Innes constant A.

        B: *float*
            Thiele-Innes constant B.

        F: *float*
            Thiele-Innes constant F.

        G: *float*
            Thiele-Innes constant G.

        periapsis_epoch: *float * or *None*
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
            self, period, eccentricity, A, B, F, G, periapsis_epoch=None,
            argument_of_latitude_reference=None, epoch_reference=None):
        self._period = period
        self._eccentricity = eccentricity
        self._Omega_node = None
        self._q = None
        self._p = None
        self._semimajor_axis = None
        self._inclination = None
        self._argument_of_latitude_reference = argument_of_latitude_reference
        self._A, self._B, self._F, self._G = A, B, F, G
        if periapsis_epoch is None:
            self._omega_periapsis = self.get_omega_periapsis()*np.pi/180.
            self._epoch_reference = epoch_reference
        self._periapsis_epoch = self._check_for_and_get_periapsis_epoch(
                periapsis_epoch, argument_of_latitude_reference, epoch_reference)

    def orbit_elements_dict_degrees(self):
        """
        Return standard orbital elements in a dictionary with angles in degrees.
        """
        dict_out = {'period': self._period, 'semimajor_axis': self.get_semimajor_axis(),
                    'eccentricity': self._eccentricity, 'inclination': self.get_inclination(),
                    'omega_periapsis': self.get_omega_periapsis(), 'Omega_node': self.get_Omega_node(),
                    'periapsis_epoch': self._periapsis_epoch
                    }
        if self._argument_of_latitude_reference is not None:
            dict_out['argument_of_latitude_reference'] = self._argument_of_latitude_reference
            dict_out['epoch_reference'] = self._epoch_reference

        return dict_out

    def thiele_innes_orbit_elements_dict_degrees(self):
        """
        Return standard orbital elements and Thiele-Innes constants in a dictionary with angles in degrees.
        """
        dict_out = self.orbit_elements_dict_degrees()
        dict_out['A'] = self._A
        dict_out['B'] = self._B
        dict_out['F'] = self._F
        dict_out['G'] = self._G
        return dict_out

    def get_semimajor_axis(self):
        """
        Return semimajor_axis in the same units as provided in input.
        """
        if self._q is None:
            self._set_q()
        if self._p is None:
            self._set_p()
        self._semimajor_axis = np.sqrt(self._p + np.sqrt(self._p**2. - self._q**2.))
        return self._semimajor_axis

    def _set_q(self):
        """
        Calculate q = A*G - B*F  = a^2*cos(i)
        """
        self._q = self._A * self._G - self._B * self._F

    def _set_p(self):
        """
        Calculate p = 0.5 * (A^2 + B^2 + F^2 + G^2) = a^2 + a^2*cos(i)^2
        """
        self._p = (self._A**2. + self._B**2. + self._F**2. + self._G**2.)/2.

    def get_inclination(self):
        """
        Return inclination in degrees. Keeping conventions inclination in range [0, 180] degrees
        """
        if self._inclination is None:
            if self._semimajor_axis is None:
                self.get_semimajor_axis()
            self._inclination = np.arccos(round(self._q / self._semimajor_axis**2., 10)) % np.pi

        return self._inclination * 180. / np.pi

    def get_omega_periapsis(self):
        """
        Return omega_periapsis in degrees. Keeping conventions Omega_node < 360 degrees
        """
        if self._Omega_node is None:
            self.get_Omega_node()
        return self._omega_periapsis * 180. / np.pi

    def get_Omega_node(self):
        """
        Return Omega_node in degrees. Keeping conventions Omega_node < 180 degrees
        """
        Omega_node = self._get_Omega_node()
        Omega_node = Omega_node % (2.*np.pi)
        if Omega_node > np.pi:
            Omega_node -= np.pi
        self._Omega_node = Omega_node
        self._omega_periapsis = self._arctan2BminusFoverAplusG - Omega_node
        self._omega_periapsis = self._omega_periapsis % (2.*np.pi)
        return self._Omega_node * 180. / np.pi

    def _get_Omega_node(self):
        """
        Calculate Omega_node from Thiele-Innes elements.
        """
        A, B, F, G = self._A, self._B, self._F, self._G
        self._arctan2BminusFoverAplusG = np.arctan2(B - F, A + G)
        return 0.5*(self._arctan2BminusFoverAplusG - np.arctan2(-B-F, A-G)) % (2.*np.pi)

    def get_reference_plane_position(self, time):
        """
        Calculate position in the reference plane at given time.
        The result is in the same units as semimajor_axis.
        """
        eccentric_anomaly = self._get_eccentric_anomaly(time)
        X, Y = self._eliptical_retangular_coordinates(eccentric_anomaly, self._eccentricity)

        matrix = np.array(([self._A, self._F],
                           [self._B, self._G]))

        projected = np.matmul(matrix, np.array([X, Y]))

        return projected

    def _eliptical_retangular_coordinates(self, eccentric_anomaly, eccentricity):
        """
        Calculate rectangular coordinates in the orbital plane
        based on eccentric anomaly and eccentricity.
        """
        X = np.cos(eccentric_anomaly) - eccentricity
        Y = np.sqrt(1 - eccentricity**2.) * np.sin(eccentric_anomaly)
        return [X, Y]
