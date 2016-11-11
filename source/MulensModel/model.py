import numpy as np
from astropy.coordinates import SkyCoord, get_body_barycentric, EarthLocation
from astropy import units as u
from astropy.time import Time, TimeDelta

from MulensModel.modelparameters import ModelParameters

class Model(object):
    """
    Caveats:
    1. Does not currently have self-consistency checks: e.g. it is
    possible to define s, a_proj, and a source distance that are not
    self-consistent. Under these circumstances, the behavior may be
    unpredictable.
    """
    def __init__(self, parameters=None,
                 t_0=0., u_0=None, t_E=0., rho=None, s=None, q=None,
                 alpha=None,
                 pi_E=None, pi_E_N=None, pi_E_E=None,
                 pi_E_ref=None,
                 lens=None, source=None, mu_rel=None,
                 coords=None, ra=None, dec=None):
        """
        Three ways to define the model:
        1. parameters = a ModelParameters() object
        2. specify t_0, u_0, t_E (optionally: rho, s, q, alpha,pi_E)
        3. specify physical properties: lens= a Lens() object, 
            source= a Source() object, mu_rel
        method 3 not implemented.

        When defining event coordinates, may specify coords as an
        astropy.coordinates.SkyCoord object, otherwise assumes RA is
        in hour angle and DEC is in degrees.
        """
        self.parameters = ModelParameters(t_0=t_0, u_0=u_0, t_E=t_E)
        if parameters is not None:
            self.parameters = ModelParameters()
        if t_0 is not None:
            self.t_0 = t_0
        if u_0 is not None:
            self.u_0 = u_0
        if t_E is not None:
            self.t_E = t_E
        if rho is not None:
            self.rho = rho
        
        par_msg = 'Must specify both or neither of pi_E_N and pi_E_E'
        if pi_E is not None:
            if pi_E_ref is None:
                self.pi_E = pi_E
            else:
                self.parameters.pi_E = MulensParallaxVector(pi_E, ref=pi_E_ref)
        if pi_E_N is not None:
            if pi_E_E is not None:
                if pi_E_ref is None:
                    self.pi_E_N = pi_E_N
                    self.pi_E_E = pi_E_E
                else:
                    self.parameters.pi_E = MulensParallaxVector(
                        pi_E_1=pi_E_N, pi_E_2=pi_E_E, ref=pi_E_ref)
            else:
                raise AttributeError(par_msg)
        else:
            if pi_E_E is not None:
                raise AttributeError(par_msg)

        if lens is not None:
            pass
        if source is not None:
            pass

        coords_msg = 'Must specify both or neither of ra and dec'
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

        self._magnification = None
        self._parallax_earth_orbital = False
        self._parallax_satellite = False
        self._parallax_topocentric = False

    @property
    def t_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self.parameters.t_0

    @t_0.setter
    def t_0(self, value):
        self.parameters.t_0 = value

    @property
    def u_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self.parameters._u_0
    
    @u_0.setter
    def u_0(self, value):
        self.parameters._u_0 = value

    @property
    def t_E(self):
        """
        The Einstein timescale. An astropy.Quantity. "day" is the default unit.
        """
        return self.parameters.t_E

    @t_E.setter
    def t_E(self, value):
        self.parameters.t_E = value

    @property
    def pi_E(self):
        """
        The microlens parallax vector. May be specified either
        relative to the sky ("NorthEast") or relative to the binary
        axis ("ParPerp"). "NorthEast" is default. A
        MulensParallaxVector object.
        """
        return self.parameters.pi_E

    @pi_E.setter
    def pi_E(self, value):
        self.parameters.pi_E = value

    @property
    def pi_E_N(self):
        """
        The North component of the microlens parallax vector.
        """
        return self.parameters.pi_E_N

    @pi_E_N.setter
    def pi_E_N(self, value):
        self.parameters.pi_E_N = value

    @property
    def pi_E_E(self):
        """
        The East component of the microlens parallax vector.
        """
        return self.parameters.pi_E_E

    @pi_E_E.setter
    def pi_E_E(self, value):
        self.parameters.pi_E_E = value

    def _trajectory(self):
        """calculates source trajectory"""
        self._trajectory_x = []
        self._trajectory_y = []
        if self._parallax_satellite is not False:
            raise ValueError('Satellite parallax not coded yet')
        if self._parallax_topocentric is not False:
            raise ValueError('Topocentric parallax not coded yet')        
        for dataset in self._datasets:
            vector_x = (dataset.time - self.t_0) / self.t_E
            vector_y = self.u_0 * np.ones(len(dataset.time))
            if self._parallax_earth_orbital is True:
                earth_center = EarthLocation.from_geocentric(0., 0., 0., u.m)
                time_ref = Time(self.t_0_par+2450000., format="jd", 
                          location=earth_center)
                dt = TimeDelta(3600.0, format='sec')
                position_ref = get_body_barycentric(body='earth', time=time_ref)
                position_ref_dt = get_body_barycentric(body='earth', time=time_ref+dt)
                velocity = (position_ref_dt.xyz - position_ref.xyz) / dt
                position = get_body_barycentric(body='earth', time=dataset._time.astropy_time)
                delta_time = dataset._time.astropy_time - time_ref
                product = np.outer(delta_time.to(u.d).value, velocity.value) * u.d * velocity.unit 
                # We calculated prodcut in this strange way because np.outer() 
                # destroys information about units of its arguments
                delta_s = position.xyz.T - product - position_ref.xyz.T
                delta_s_e = -delta_s[:,0]
                delta_s_n = delta_s[:,2]
                delta_tau = delta_s_n.value * self.pi_E_N + delta_s_e.value * self.pi_E_E
                delta_beta = -delta_s_n.value * self.pi_E_E + delta_s_e.value * self.pi_E_N
                vector_x += delta_tau
                vector_y += delta_beta
            self._trajectory_x.append(vector_x)
            self._trajectory_y.append(vector_y)

    @property
    def magnification(self):
        """a list of magnifications calculated for every dataset time vector"""
        if self._magnification is not None:
            return self._magnification
        self._magnification = []
        self._trajectory()
        for i_data in range(len(self._datasets)):
            u2 = self._trajectory_x[i_data]**2 + self._trajectory_y[i_data]**2
            self._magnification.append((u2 + 2.) / np.sqrt(u2 * (u2 + 4.)))
        return self._magnification

    @magnification.setter
    def magnification(self, new_value):
        self._magnification = new_value

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
        """
        Does not work. See __name__=="__main__"
        """
        try:
            self._coords.ra = new_value
        except AttributeError:
            self._coords = SkyCoord(new_value, 0.0, unit=(u.hourangle, u.deg))

    @property
    def dec(self):
        """
        Declination
        """
        return self._coords.dec

    @dec.setter
    def dec(self, new_value):
        """
        Does not work. See __name__=="__main__"
        """
        try:
            self._coords.dec = new_value
        except AttributeError:
            self._coords = SkyCoord(0.0, new_value, unit=(u.hourangle, u.deg))

    def parallax(self, earth_orbital=False, satellite=False, topocentric=False):
        """specifies which types of the parallax will be included in calculations"""
        self._parallax_earth_orbital = earth_orbital
        self._parallax_satellite = satellite
        self._parallax_topocentric = topocentric

if __name__ == "__main__":
    model_1 = Model(coords="18:00:00 -30:00:00")
    print(model_1.coords,model_1.ra,model_1.dec)

    model_2 = Model()
    model_2.ra = "17:00:00"
    print(model_2.coords,model_2.ra,model_2.dec)
    model_2.dec = "40:03:01"
    print(model_2.coords,model_2.ra,model_2.dec)

    model_3 = Model()
    model_3.coords = "17:00:00 -27:32:14"
    print(model_3.coords,model_3.ra,model_3.dec)

