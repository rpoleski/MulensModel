import numpy as np
from math import fsum
from astropy.coordinates import SkyCoord, get_body_barycentric, EarthLocation
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy.coordinates import GeocentricTrueEcliptic
from astropy import units as u
from astropy import _erfa as erfa

from MulensModel.modelparameters import ModelParameters
from MulensModel.mulenstime import MulensTime


def dot(cartesian, vector):
    """dot product of Astropy CartersianRepresentation and np.array"""
    return cartesian.x * vector[0] + cartesian.y * vector[1] + cartesian.z * vector[2]


class Model(object):
    """
    Caveats:
    1. Does not currently have self-consistency checks: e.g. it is
    possible to define s, a_proj, and a source distance that are not
    self-consistent. Under these circumstances, the behavior may be
    unpredictable.
    """
    def __init__(self, parameters=None,
                 t_0=None, u_0=None, t_E=None, rho=None, s=None, q=None,
                 alpha=None,
                 pi_E=None, pi_E_N=None, pi_E_E=None,
                 pi_E_ref=None, t_0_par=None, 
                 lens=None, source=None, mu_rel=None,
                 coords=None, ra=None, dec=None):
        """
        Three ways to define the model:
        1. parameters = a ModelParameters() object
        2. specify t_0, u_0, t_E (optionally: rho, s, q, alpha, pi_E, t_0_par)
        3. specify physical properties: lens= a Lens() object, 
            source= a Source() object, mu_rel
        method 3 not implemented.

        When defining event coordinates, may specify coords as an
        astropy.coordinates.SkyCoord object, otherwise assumes RA is
        in hour angle and DEC is in degrees.
        """
        if isinstance(parameters, ModelParameters):
            self._parameters = parameters
        elif parameters is None:
            self._parameters = ModelParameters()
        else:
            raise TypeError(
                "If specified, parameters must be a ModelParameters object.")
        if t_0 is not None:
            self.t_0 = t_0
        if u_0 is not None:
            self.u_0 = u_0
        if t_E is not None:
            self.t_E = t_E
        if rho is not None:
            self.rho = rho
        self.t_0_par = t_0_par
        
        par_msg = 'Must specify both or neither of pi_E_N and pi_E_E'
        if pi_E is not None:
            if pi_E_ref is None:
                self.pi_E = pi_E
            else:
                self._parameters.pi_E = MulensParallaxVector(pi_E, ref=pi_E_ref)
        if pi_E_N is not None:
            if pi_E_E is not None:
                if pi_E_ref is None:
                    self.pi_E_N = pi_E_N
                    self.pi_E_E = pi_E_E
                else:
                    self._parameters.pi_E = MulensParallaxVector(
                        pi_E_1=pi_E_N, pi_E_2=pi_E_E, ref=pi_E_ref)
            else:
                raise AttributeError(par_msg)
        else:
            if pi_E_E is not None:
                raise AttributeError(par_msg)

        #Replace with MulensSystem object
        #If lens and source are defined, but t0, tE are not, then 
        #set t0=0.,tE=1.
        if lens is not None :
            if source is None:
                raise AttributeError(
                    'If lens is specified, source must also be specified.')
            else:
                self.lens = lens
                self.source = source
        if source is not None and lens is None:
            raise AttributeError(
                'If source is specified, lens must also be specified.')

        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._coords = coords
            else:
                self._coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        if ra is not None:
            if dec is not None:
                self._coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            else:
                raise AttributeError(coords_msg)
        else:
            if ra is not None:
                raise AttributeError(coords_msg)

        self.reset_magnification()
        self._parallax_earth_orbital = False
        self._parallax_satellite = False
        self._parallax_topocentric = False
        self._delta_annual = {}
        self._delta_satellite = {}

    @property
    def t_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self._parameters.t_0

    @t_0.setter
    def t_0(self, value):
        self._parameters.t_0 = value
        self.reset_magnification()

    @property
    def u_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self._parameters._u_0
    
    @u_0.setter
    def u_0(self, value):
        self._parameters._u_0 = value
        self.reset_magnification()

    @property
    def t_E(self):
        """
        The Einstein timescale. An astropy.Quantity. "day" is the default unit.
        """
        return self._parameters.t_E

    @t_E.setter
    def t_E(self, value):
        self._parameters.t_E = value
        self.reset_magnification()
        
    @property
    def pi_E(self):
        """
        The microlens parallax vector. May be specified either
        relative to the sky ("NorthEast") or relative to the binary
        axis ("ParPerp"). "NorthEast" is default. A
        MulensParallaxVector object.
        """
        return self._parameters.pi_E

    @pi_E.setter
    def pi_E(self, value):
        self._parameters.pi_E = value
        self.reset_magnification()

    @property
    def pi_E_N(self):
        """
        The North component of the microlens parallax vector.
        """
        return self._parameters.pi_E_N

    @pi_E_N.setter
    def pi_E_N(self, value):
        self._parameters.pi_E_N = value
        self.reset_magnification()

    @property
    def pi_E_E(self):
        """
        The East component of the microlens parallax vector.
        """
        return self._parameters.pi_E_E

    @pi_E_E.setter
    def pi_E_E(self, value):
        self._parameters.pi_E_E = value
        self.reset_magnification()

    @property
    def t_0_par(self):
        """reference time for parameters, in particular microlensing parallax"""
        return self._t_0_par

    @t_0_par.setter
    def t_0_par(self, value):
        if isinstance(value, MulensTime):
            self._t_0_par = value
        elif value is None:
            self._t_0_par = None
        else:
            self._t_0_par = MulensTime(value)
        self.reset_magnification()

    def parameters(
        self, t_0=0., u_0=None, t_E=1., rho=None, s=None, q=None, alpha=None, 
        pi_E=None, pi_E_N=None, pi_E_E=None, pi_E_ref=None):
        if u_0 is None:
            return "{0}".format(self._parameters)
        else:
            self._parameters = ModelParameters(
                t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, s=s, q=q, alpha=alpha, 
                pi_E=pi_E, pi_E_N=pi_E_N, pi_E_E=pi_E_E, pi_E_ref=pi_E_ref)

    def _trajectory(self):
        """calculates source trajectory"""
        self._trajectory_x = []
        self._trajectory_y = []
        if self._parallax_topocentric is not False:
            raise ValueError('Topocentric parallax not coded yet')
        n_satellite = 0 
        for dataset in self._datasets:
            vector_tau = (dataset.time + dataset.time_zeropoint - (self.t_0 + self._parameters._t_0.zeropoint)) / self.t_E
            vector_u = self.u_0 * np.ones(dataset.n_epochs)
            if self._parallax_earth_orbital:
                (delta_tau, delta_u) = self._annual_parallax_trajectory(dataset)
                vector_tau += delta_tau
                vector_u += delta_u
            if self._parallax_satellite and dataset.is_satellite: 
                (delta_tau, delta_u) = self._satellite_parallax_trajectory(dataset)
                vector_tau += delta_tau
                vector_u += delta_u
                n_satellite += 1
            if self._parameters.n_lenses == 1:
                vector_x = vector_tau
                vector_y = vector_u
            elif self._parameters.n_lenses == 2:
                sin_alpha = np.sin(self._parameters.alpha)
                cos_alpha = np.cos(self._parameters.alpha)
                vector_x = vector_u * sin_alpha - vector_tau * cos_alpha
                vector_y = -vector_u * cos_alpha - vector_tau * sin_alpha
                vector_x += self._parameters._s / 2.
            else:
                raise Exception("trajectory for more than 2 lenses not handled yet")
            self._trajectory_x.append(vector_x)
            self._trajectory_y.append(vector_y)
        if self._parallax_satellite and n_satellite == 0:
            raise ValueError('Satellite parallax turned on, but no satellite data provided')

    def _get_delta_satellite(self, dataset):
        """calculates differences of Earth and satellite positions projected on the plane of the sky at event position"""
        direction = np.array(self._coords.cartesian.xyz.value)
        north = np.array([0., 0., 1.])
        east_projected = np.cross(north, direction)
        east_projected /= np.linalg.norm(east_projected)
        north_projected = np.cross(direction, east_projected)
        satellite = dataset.satellite_skycoord
        satellite.transform_to(frame=self._coords.frame) # We want to be sure frames are the same.
        self._delta_satellite[dataset] = {}
        self._delta_satellite[dataset]['N'] = -dot(satellite.cartesian, north_projected)
        self._delta_satellite[dataset]['E'] = -dot(satellite.cartesian, east_projected)
        self._delta_satellite[dataset]['D'] = -dot(satellite.cartesian, direction)

    def _satellite_parallax_trajectory(self, dataset):
        """calcualate satellite parallax component of trajectory"""
        if dataset not in self._delta_satellite:
            self._get_delta_satellite(dataset)
        delta_tau = self._delta_satellite[dataset]['N'] * self.pi_E_N + self._delta_satellite[dataset]['E'] * self.pi_E_E
        delta_beta = -self._delta_satellite[dataset]['N'] * self.pi_E_E + self._delta_satellite[dataset]['E'] * self.pi_E_N
        return (delta_tau, delta_beta)

    def _get_delta_annual(self, dataset):
        """calculates projected Earth positions required by annual parallax"""
        if self.t_0_par is None:
            msg1 = 'Annual parallax effect cannot be '
            msg2 = 'calculated if t_0_par is not set'
            raise ValueError(msg1 + msg2)
        time_ref = self.t_0_par.astropy_time
        position_ref = get_body_barycentric(body='earth', time=time_ref)
        # the 3 lines below, that calculate velocity for t_0_par, are based on astropy 1.3 
        # https://github.com/astropy/astropy/blob/master/astropy/coordinates/solar_system.py
        jd1, jd2 = get_jd12(time_ref, 'tdb')
        earth_pv_helio, earth_pv_bary = erfa.epv00(jd1, jd2)
        velocity = earth_pv_bary[..., 1, :] * u.au / u.day
        
        position = get_body_barycentric(body='earth', time=dataset._time.astropy_time)
        delta_time = dataset._time.astropy_time - time_ref
        product = np.outer(delta_time.to(u.d).value, velocity.value) * u.d * velocity.unit
        # We calculated product in this strange way because np.outer()
        # destroys information about units of its arguments.
        delta_s = position.xyz.T - product - position_ref.xyz.T
        self._delta_annual[dataset] = {}
        self._delta_annual[dataset]['E'] = -delta_s[:,0]
        self._delta_annual[dataset]['N'] = delta_s[:,2]

    def _annual_parallax_trajectory(self, dataset):
        """calcualate annual parallax component of trajectory"""
        if dataset not in self._delta_annual:
            self._get_delta_annual(dataset)
        delta_tau = self._delta_annual[dataset]['N'] * self.pi_E_N + self._delta_annual[dataset]['E'] * self.pi_E_E
        delta_beta = -self._delta_annual[dataset]['N'] * self.pi_E_E + self._delta_annual[dataset]['E'] * self.pi_E_N
        return (delta_tau, delta_beta)

    @property
    def magnification(self):
        """a list of magnifications calculated for every dataset time vector"""
        if self._magnification is not None:
            return self._magnification
        self._magnification = []
        self._trajectory()
        for i_data in range(len(self._datasets)):
            if self._parameters.n_lenses == 1:
                u2 = self._trajectory_x[i_data]**2 + self._trajectory_y[i_data]**2
                self._magnification.append((u2 + 2.) / np.sqrt(u2 * (u2 + 4.)))
            elif self._parameters.n_lenses == 2:
                q = self._parameters.q
                s = self._parameters.s
                m1 = 1. / (1. + q)
                m2 = q / (1. + q)
                # See Witt & Mao (1995) for details.
                m = 0.5 * (m1 + m2) # Quantity m equals 1/2 according to the above definition
                dm = 0.5 * (m2 - m1)
                z1 = -0.5 * s + 0.j
                z2 = 0.5 * s + 0.j
                zeta = self._trajectory_x[i_data] - s + self._trajectory_y[i_data] * 1.j
                zeta_c = np.conjugate(zeta)
                c5 = z1**2 - zeta_c**2
                c4 = -2.*m*zeta_c + zeta*zeta_c**2 - 2.*dm*z1 - zeta*z1**2
                c3 = 4.*m*zeta*zeta_c + 4.*dm*zeta_c*z1 + 2.*(zeta_c*z1)**2 - 2.*z1**4
                c2 = 4.*m**2*zeta + 4.*m*dm*z1 - 4.*dm*zeta*zeta_c*z1 - 2.*zeta*zeta_c**2*z1**2 + 4.*dm*z1**3 + 2.*zeta*z1**4
                c1 = -8.*m*dm*zeta*z1 - 4.*dm**2*z1**2 - 4.*(m*z1)**2 - 4.*m*zeta*zeta_c*z1**2 - 4.*dm*zeta_c*z1**3 - zeta_c**2*z1**4 + z1**6
                c0 = 4.*dm**2*zeta + 4.*m*dm*z1 + 4.*dm*zeta*zeta_c*z1 + 2.*m*zeta_c*z1**2 + zeta*zeta_c**2*z1**2 - 2.*dm*z1**3 - zeta*z1**4
                c0 *= z1**2
                c = np.array([c0, c1, c2, c3, c4, c5]).reshape(6)
                roots = np.polynomial.polynomial.polyroots(c)
                solution_c = np.conjugate(zeta_c + m1 / (roots-z1) + m2 / (roots-z2))
                correct = np.ones(5, dtype = bool)
                for i in range(5): # Here we check if solutions of the polynomial are also solutions of the lens equation
                    j = np.argmin(abs((solution_c-roots[i])**2))
                    print(i, abs((solution_c-roots[i])**2)[j], roots[i], solution_c[j])
                    if i != j:
                        correct[i] = False
                roots_c = np.conjugate(roots)
                dzeta_dz_c = m1 / (z1-roots_c)**2 + m2 / (z2-roots_c)**2
                det_j = 1. - dzeta_dz_c * np.conjugate(dzeta_dz_c)
                self._magnification.append([fsum(1./abs(det_j[correct]))])
            else:
                raise Exception("magnification for more than 2 lenses not handled yet")
        return self._magnification

    @magnification.setter
    def magnification(self, new_value):
        self._magnification = new_value
        
    def reset_magnification(self):
        """destroy existing information on magnification - call it after you change parameters"""
        self._magnification = None

    def set_datasets(self, datasets):
        """set _datasets property"""
        self._datasets = datasets

    @property
    def lens(self):
        pass

    @lens.setter
    def lens(self, new_lens):
        pass

    @property
    def source(self):
        pass

    @source.setter
    def source(self, new_lens):
        pass

    @property
    def coords(self):
        """
        Sky coordinates (RA,Dec)
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        if isinstance(new_value, SkyCoord):
            self._coords = new_value
        else:
            self._coords = SkyCoord(new_value, unit=(u.hourangle, u.deg))

    @property
    def ra(self):
        """
        Right Ascension
        """
        return self._coords.ra

    @ra.setter
    def ra(self, new_value):
        try:
            self._coords.ra = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    new_value, 0.0, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    new_value, self._coords.dec, unit=(u.hourangle, u.deg)) 

    @property
    def dec(self):
        """
        Declination
        """
        return self._coords.dec

    @dec.setter
    def dec(self, new_value):
        try:
            self._coords.dec = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    0.0, new_value, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    self._coords.ra, new_value, unit=(u.hourangle, u.deg))

    @property
    def galactic_l(self):
        """Galactic longitude"""
        l = self._coords.galactic.l
        if l > 180.*u.deg:
            l = l - 360*u.deg
        return l

    @property
    def galactic_b(self):
        """Galactic latitude"""
        return self._coords.galactic.b

    @property
    def ecliptic_lon(self):
        """ecliptic longitude"""
        return self._coords.transform_to(GeocentricTrueEcliptic).lon

    @property
    def ecliptic_lat(self):
        """ecliptic latitude"""
        return self._coords.transform_to(GeocentricTrueEcliptic).lat

    def parallax(self, earth_orbital=False, satellite=False, topocentric=False):
        """specifies which types of the parallax will be included in calculations"""
        self._parallax_earth_orbital = earth_orbital
        self._parallax_satellite = satellite
        self._parallax_topocentric = topocentric

