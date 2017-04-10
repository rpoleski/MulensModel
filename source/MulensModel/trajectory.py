import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric, EarthLocation
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
           parallax - boolean dictionary specifying what parallax effects should be used.
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
        self.parallax = parallax
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
        
        #Apply Earth Orbital parallax effect
        if self.parallax['earth_orbital']:
            [delta_tau, delta_u] = self._annual_parallax_trajectory()
            vector_tau += delta_tau
            vector_u += delta_u

        #Apply satellite parallax effect
        if self.parallax['satellite'] and self.satellite_skycoord is not None: 
            [delta_tau, delta_u] = self._satellite_parallax_trajectory()
            vector_tau += delta_tau
            vector_u += delta_u

        #If 2 lenses, rotate trajectory relative to binary axis
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
                      + delta['E'] * self.parameters. pi_E_E)
        delta_beta = ( -delta['N'] * self.parameters.pi_E_E 
                      + delta['E'] * self.parameters.pi_E_N)
        return [delta_tau, delta_beta]

    def _annual_parallax_trajectory(self):
        """calcualate annual parallax component of trajectory"""
        print('WARNING - Some probable problems with annual parallax due to interface with astropy functions get_body_barycentric and get_jd12. These must depend on both reference frame and time standard. But it is only possible to set time standard and jd vs. mjd. - JCY')

        #Calculate the parallax offsets
        delta_annual = self._get_delta_annual()
        return self._project_delta(delta_annual)

    def _get_delta_annual(self):
        """
        calculates projected Earth positions required by annual parallax
        """
        print('WARNING: Does not take into account coordinates!!!!')

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
        jd1, jd2 = get_jd12(Time(time_ref,format='jd',scale='tdb'), 'tdb')
        earth_pv_helio, earth_pv_bary = erfa.epv00(jd1, jd2)
        velocity = earth_pv_bary[..., 1, :] * u.au / u.day
        
        position = get_body_barycentric(
            body='earth', time=Time(self.times,format='jd', scale='tdb'))
        delta_time = self.times - time_ref
        product = (np.outer(delta_time, velocity.value) * u.d * velocity.unit)
        # We calculated product in this strange way because np.outer()
        # destroys information about units of its arguments.
        delta_s = position.xyz.T - product - position_ref.xyz.T

        #Return annual parallax offsets
        return {'N': delta_s[:,2], 'E': -delta_s[:,0]}


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
        delta_satellite['N'] = -dot(satellite.cartesian, north_projected)
        delta_satellite['E'] = -dot(satellite.cartesian, east_projected)
        delta_satellite['D'] = -dot(satellite.cartesian, direction)

        return delta_satellite

