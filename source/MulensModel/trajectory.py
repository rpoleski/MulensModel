import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric, EarthLocation
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy import _erfa as erfa
from astropy.time import Time

class Trajectory(object):
    """
    The (dimensionless) X, Y trajectory of the source in the
    source plane. t_0_par is the reference time for the parallax
    vector. If not set, defaults to t_0.
    """
    def __init__(
        self, times, parameters=None, parallax=None, t_0_par=None,
        coords=None, satellite_coords=None):

        #Save parameters
        if isinstance(times, (list, tuple, np.ndarray)):
            self.times = times
        else:
            self.times = np.array(times)
        self.parameters = parameters
        self.parallax = parallax
        self.t_0_par = t_0_par
        self.coords = coords
        self.satellite_coords = satellite_coords

        #Calculate trajectory
        self.get_xy()

    def get_xy(self):
        """
        For a given set of parameters (a ModelParameters object),
        calculate the xy position of the source.
        """
        vector_tau = (
            (self.times - self.parameters.t_0)
            / float(self.parameters.t_E))
        vector_u = self.parameters.u_0 * np.ones(self.times.size)
        
        if self.parallax['earth_orbital']:
            [delta_tau, delta_u] = self._annual_parallax_trajectory()
            vector_tau += delta_tau
            vector_u += delta_u


        if self.parallax['satellite'] and self.satellite_coords is not None: 
            [delta_tau, delta_u] = self._satellite_parallax_trajectory()
            vector_tau += delta_tau
            vector_u += delta_u
            n_satellite += 1

        if self.parameters.n_lenses == 1:
            vector_x = vector_tau
            vector_y = vector_u
        elif self.parameters.n_lenses == 2:
            sin_alpha = np.sin(self.parameters.alpha)
            cos_alpha = np.cos(self.parameters.alpha)
            vector_x = vector_u * sin_alpha - vector_tau * cos_alpha
            vector_y = -vector_u * cos_alpha - vector_tau * sin_alpha
            vector_x += self.parameters.s / 2.
        else:
            raise Exception(
                "trajectory for more than 2 lenses not handled yet")

        self.x = vector_x
        self.y = vector_y

    def _project_delta(self, delta):
        delta_tau =  (  delta['N'] * self.parameters.pi_E_N 
                      + delta['E'] * self.parameters. pi_E_E)
        delta_beta = ( -delta['N'] * self.parameters.pi_E_E 
                      + delta['E'] * self.parameters.pi_E_N)
        return [delta_tau, delta_beta]

    def _annual_parallax_trajectory(self):
        """calcualate annual parallax component of trajectory"""
        print('WARNING - Some probable problems with annual parallax due to interface with astropy functions get_body_barycentric and get_jd12. These must depend on both reference frame and time standard. But it is only possible to set time standard and jd vs. mjd. - JCY')

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

        return {'N': delta_s[:,2], 'E': -delta_s[:,0]}


    def _satellite_parallax_trajectory(self):
        """calcualate satellite parallax component of trajectory"""
        
        delta_satellite = self._get_delta_satellite()
        return self._project_delta(delta_satellite)

        
    def _get_delta_satellite(self):
        """
        calculates differences of Earth and satellite positions
        projected on the plane of the sky at event position
        """
        direction = np.array(self.coords.cartesian.xyz.value)
        north = np.array([0., 0., 1.])
        east_projected = np.cross(north, direction)
        east_projected /= np.linalg.norm(east_projected)
        north_projected = np.cross(direction, east_projected)
        satellite = self.satellite_skycoord
        # We want to be sure frames are the same.
        satellite.transform_to(frame=self._coords.frame) 
            
        delta_satellite = {}
        delta_satellite['N'] = -dot(satellite.cartesian, north_projected)
        delta_satellite['E'] = -dot(satellite.cartesian, east_projected)
        delta_satellite['D'] = -dot(satellite.cartesian, direction)

        return delta_satellite

