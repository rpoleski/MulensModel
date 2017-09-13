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
    source plane. t_0_par is the reference time for the parallax
    vector. If not set, defaults to t_0.
    For binary lens, the origin of the coordinate system is at 
    the center of mass with higher mass at negative X and Y=0.
    """
    def __init__(
        self, times, parameters=None, parallax=None, t_0_par=None,
        coords=None, satellite_skycoord=None):
        """
        Required arguments: 
           times - the times at which to generate the source trajectory, e.g. a vector.
           parameters - a ModelParameters object specifying the microlensing parameters

        Optional parallax keywords:
           parallax - boolean dictionary specifying what parallax effects should be used. Default is False. (differs from Modely.py which defaults to True)
           t_0_par - reference time for the parallax
           coords - sky coordinates of the event
           satellite_skycoord - sky coordinates of the satellite
               specified by the ephemrides file. see
               MulensData.satellite_skycoord.
        """
        #Set times
        if isinstance(times, (list, tuple, np.ndarray)):
            self.times = times
        else:
            self.times = np.array(times)

        #Check for ModelParameters and set.
        if isinstance(parameters, ModelParameters):
            self.parameters = parameters
        else:
            raise ValueError('parameters is a required and must be a ModelParameters object')

        #Set parallax values
        self.parallax = {'earth_orbital':False, 
                         'satellite':False, 
                         'topocentric':False}
        if parallax is not None:
            for (key, value) in parallax.items():
                self.parallax[key] = value

        self.t_0_par = t_0_par
        self.coords = coords
        self.satellite_skycoord = satellite_skycoord

        #Calculate trajectory
        self.get_xy()

    def get_xy(self):
        """
        For a given set of parameters (a ModelParameters object),
        calculate the xy position of the source.
        """
        #Calculate the position of the source
        vector_tau = (
            (self.times - self.parameters.t_0)
            / float(self.parameters.t_E))
        vector_u = self.parameters.u_0 * np.ones(self.times.size)
        
        #If parallax is non-zero, apply parallax effects:
        if self.parameters.pi_E is not None:
            #Apply Earth Orbital parallax effect
            if self.parallax['earth_orbital']:
                [delta_tau, delta_u] = self._annual_parallax_trajectory()
                vector_tau += delta_tau
                vector_u += delta_u

            #Apply satellite parallax effect
            if (self.parallax['satellite'] 
                and self.satellite_skycoord is not None): 
                [delta_tau, delta_u] = self._satellite_parallax_trajectory()
                vector_tau += delta_tau
                vector_u += delta_u

        #If 2 lenses, rotate trajectory relative to binary lens axis
        if self.parameters.n_lenses == 1:
            vector_x = vector_tau
            vector_y = vector_u
        elif self.parameters.n_lenses == 2:
            sin_alpha = np.sin(self.parameters.alpha)
            cos_alpha = np.cos(self.parameters.alpha)
            shift_x = - self.parameters.s * self.parameters.q / (1. + self.parameters.q)
            vector_x = vector_u * sin_alpha - vector_tau * cos_alpha + shift_x
            vector_y = -vector_u * cos_alpha - vector_tau * sin_alpha
        else:
            raise Exception(
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
        if self.t_0_par is not None:
            time_ref = self.t_0_par
        else:
            time_ref = self.parameters.t_0

        position_ref = get_body_barycentric(
            body='earth', time=Time(time_ref,format='jd',scale='tdb')) 
        #seems that get_body_barycentric depends on time system, but there is
        #no way to set BJD_TDB in astropy.Time()
        #Likewise, get_jd12 depends on time system. 

        """
        the 3 lines below, that calculate velocity for t_0_par, are
        based on astropy 1.3
        https://github.com/astropy/astropy/blob/master/astropy/coordinates/solar_system.py
        """
        (jd1, jd2) = get_jd12(Time(time_ref,format='jd',scale='tdb'), 'tdb')
        (earth_pv_helio, earth_pv_bary) = erfa.epv00(jd1, jd2)
        velocity = earth_pv_bary[..., 1, :] # This is in (u.au/u.day) 
        # but we don't multiply by unit here, because np.outer() (used later)
        # destroys information of argument units.
        
        position = get_body_barycentric(
            body='earth', time=Time(self.times,format='jd', scale='tdb'))
        product = np.outer(self.times - time_ref, velocity) * u.au
        delta_s = position.xyz.T - product - position_ref.xyz.T

        north = np.array([0., 0., 1.])
        direction = np.array(self.coords.cartesian.xyz.value)
        east_projected = utils.Utils.vector_product_normalized(north, direction)
        north_projected = utils.Utils.vector_product_normalized(direction, east_projected)
        out_n = -np.dot(delta_s.value, north_projected)
        out_e = -np.dot(delta_s.value, east_projected)
        
        return {'N': out_n, 'E': out_e}

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

        return delta_satellite
