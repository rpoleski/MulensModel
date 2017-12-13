import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric
#from astropy.coordinates import SkyCoord, EarthLocation
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy import _erfa as erfa
from astropy.time import Time

from MulensModel import utils
from MulensModel.modelparameters import ModelParameters


class Trajectory(object):
    """
    The (dimensionless) X, Y trajectory of the source in the
    source plane. This class includes internal functions that calculate 
    how microlensing parallax affects the trajectory.

    For binary lens, the origin of the coordinate system is at 
    the center of mass with higher mass at negative X and Y=0.

    Arguments : 
        times: [*float*, *list*, *np.ndarray*], required 
            the times at which to generate the source trajectory, e.g. a vector.

        parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`, required
            a ModelParameters object specifying the microlensing parameters

        parallax: *boolean dictionary*, optional
            specifies what parallax effects should be used. Default is
            *False* for each of *'earth_orbital'*, *'satellite'*, and 
            *'topocentric'*. (differs from 
            :py:class:`~MulensModel.model.Model` which defaults to *True*)

        coords: :py:class:`~MulensModel.coordinates.Coordinates`, optional
            sky coordinates of the event; required for parallax calculations 

        satellite_skycoord: *Astropy.coordinates.SkyCord*, optional 
            sky coordinates of the satellite specified by the
            ephemrides file. see
            :py:obj:`MulensModel.mulensdata.MulensData.satellite_skycoord.`
    """

    _get_delta_annual_results = dict()
    _get_delta_satellite_results = dict()

    def __init__(self, times, parameters, parallax=None,
                coords=None, satellite_skycoord=None, earth_coords=None):
        #Set times
        if isinstance(times, (list, tuple, np.ndarray)):
            self.times = times
        else:
            self.times = np.array(times)

        #Check for ModelParameters and set.
        if isinstance(parameters, ModelParameters):
            self.parameters = parameters
        else:
            m = 'parameters is a required and must be a ModelParameters object'
            raise TypeError(m)

        #Set parallax values
        self.parallax = {'earth_orbital': False, 
                         'satellite': False, 
                         'topocentric': False}
        if parallax is not None:
            for (key, value) in parallax.items():
                self.parallax[key] = value

        self.coords = coords
        self.satellite_skycoord = satellite_skycoord
        if earth_coords is not None:
            raise NotImplementedError("The earth_coords needed for " +
                    "topocentric parallax is not implemented yet")
        self.earth_coords = None

        # Calculate trajectory
        self.get_xy()

    def get_xy(self):
        """
        For a given set of parameters (a
        :py:class:`~MulensModel.modelparameters.ModelParameters`
        object), calculate the xy position of the source.
        """
        # Calculate the position of the source
        vector_tau = ((self.times - self.parameters.t_0)
                                                / float(self.parameters.t_E))
        vector_u = self.parameters.u_0 * np.ones(self.times.size)
        
        # If parallax is non-zero, apply parallax effects:
        if self.parameters.pi_E is not None:
            # Apply Earth Orbital parallax effect
            if self.parallax['earth_orbital']:
                [delta_tau, delta_u] = self._annual_parallax_trajectory()
                vector_tau += delta_tau
                vector_u += delta_u

            # Apply satellite parallax effect
            if (self.parallax['satellite'] 
                                    and self.satellite_skycoord is not None): 
                [delta_tau, delta_u] = self._satellite_parallax_trajectory()
                vector_tau += delta_tau
                vector_u += delta_u

            # Apply topocentric parallax effect
            if self.parallax['topocentric'] and self.earth_coords is not None:
                # When you implemenet it, make sure the behaviour depends on 
                # the access to the observatory location information as 
                # the satellite parallax depends on the acces to 
                # satellite_skycoord.
                raise NotImplementedError("The topocentric parallax effect " +
                                          "not implemented yet")

        # If 2 lenses, rotate trajectory relative to binary lens axis
        if self.parameters.n_lenses == 1:
            vector_x = vector_tau
            vector_y = vector_u
        elif self.parameters.n_lenses == 2:
            sin_alpha = np.sin(np.deg2rad(self.parameters.alpha))
            cos_alpha = np.cos(np.deg2rad(self.parameters.alpha))
            shift_x = - self.parameters.s * self.parameters.q / (1. +
                                                            self.parameters.q)
            vector_x = vector_u * sin_alpha - vector_tau * cos_alpha + shift_x
            vector_y = -vector_u * cos_alpha - vector_tau * sin_alpha
            # The above equations use alpha in counterclockwise convention, 
            # i.e., the same as proposed by Skowron et al. (2011)
        else:
            raise NotImplementedError(
                "trajectory for more than 2 lenses not handled yet")

        #Store trajectory
        self.x = vector_x
        self.y = vector_y

    def _project_delta(self, delta):
        """
        Project N and E parallax offset vector onto the tau, beta plane.
        """
        delta_tau =  (  delta['N'] * self.parameters.pi_E_N 
                      + delta['E'] * self.parameters.pi_E_E)
        delta_beta = ( -delta['N'] * self.parameters.pi_E_E 
                      + delta['E'] * self.parameters.pi_E_N)
        return [delta_tau, delta_beta]

    def _annual_parallax_trajectory(self):
        """calcualate annual parallax component of trajectory"""

        #Calculate the parallax offsets
        delta_annual = self._get_delta_annual()
        return self._project_delta(delta_annual)

    def _get_delta_annual(self):
        """
        calculates projected Earth positions required by annual parallax
        """
        index = (self.parameters.t_0_par, hash(self.coords), tuple(self.times.tolist()))
        if index in self._get_delta_annual_results:
            return self._get_delta_annual_results[index]
        time_ref = self.parameters.t_0_par

        position_ref = get_body_barycentric(
            body='earth', time=Time(time_ref, format='jd',scale='tdb')) 
        #seems that get_body_barycentric depends on time system, but there is
        #no way to set BJD_TDB in astropy.Time()
        #Likewise, get_jd12 depends on time system. 

        """
        the 3 lines below, that calculate velocity for t_0_par, are
        based on astropy 1.3
        https://github.com/astropy/astropy/blob/master/astropy/coordinates/solar_system.py
        """
        (jd1, jd2) = get_jd12(Time(time_ref, format='jd', scale='tdb'), 'tdb')
        (earth_pv_helio, earth_pv_bary) = erfa.epv00(jd1, jd2)
        velocity = earth_pv_bary[..., 1, :] # This is in (u.au/u.day) 
        # but we don't multiply by unit here, because np.outer() (used later)
        # destroys information of argument units.
        
        if not np.all(np.isfinite(self.times)):
            msg = "Some times have incorrect values: {:}".format(
                    self.times[~np.isfinite(self.times)])
            raise ValueError(msg)

        position = get_body_barycentric(
            body='earth', time=Time(self.times, format='jd', scale='tdb'))
        product = np.outer(self.times - time_ref, velocity) * u.au
        delta_s = position.xyz.T - product - position_ref.xyz.T

        north = np.array([0., 0., 1.])
        direction = np.array(self.coords.cartesian.xyz.value)
        vector_product_normalized = utils.Utils.vector_product_normalized
        east_projected = vector_product_normalized(north, direction)
        north_projected = vector_product_normalized(direction, east_projected)
        out_n = -np.dot(delta_s.value, north_projected)
        out_e = -np.dot(delta_s.value, east_projected)

        out = {'N': out_n, 'E': out_e}
        self._get_delta_annual_results[index] = out
        return out

    def _satellite_parallax_trajectory(self):
        """calcualate satellite parallax component of trajectory"""
        
        #Calculate the parallax offset due to the satellite
        delta_satellite = self._get_delta_satellite()
        return self._project_delta(delta_satellite)
 
    def _get_delta_satellite(self):
        """
        calculates differences of Earth and satellite positions
        projected on the plane of the sky at event position
        """
        index = (hash(self.coords), hash(self.satellite_skycoord))
        if index in self._get_delta_satellite_results.keys():
            return self._get_delta_satellite_results[index]

        #Set the N,E coordinate frame based on the direction of the event
        direction = np.array(self.coords.cartesian.xyz.value)
        north = np.array([0., 0., 1.])
        east_projected = np.cross(north, direction)
        east_projected /= np.linalg.norm(east_projected)
        north_projected = np.cross(direction, east_projected)

        #Transform the satellite ephemrides to that coordinate system
        satellite = self.satellite_skycoord 
        satellite.transform_to(frame=self.coords.frame) 

        #Project the satellite parallax effect
        delta_satellite = {}
        dot = utils.Utils.dot
        delta_satellite['N'] = -dot(satellite.cartesian, north_projected).value
        delta_satellite['E'] = -dot(satellite.cartesian, east_projected).value
        delta_satellite['D'] = -dot(satellite.cartesian, direction).value

        self._get_delta_satellite_results[index] = delta_satellite
        return delta_satellite

