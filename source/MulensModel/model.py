import numpy as np
from math import fsum
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import GeocentricTrueEcliptic
import matplotlib.pyplot as pl

from MulensModel.modelparameters import ModelParameters
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData

#JCY: some probable problems with annual parallax due to interface
#with astropy functions get_body_barycentric and get_jd12. These must
#depend on both reference frame and time standard. But it is only
#possible to set time standard and jd vs. mjd.

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

    2. satellite parallax works for datasets, but not for
    model. i.e. The satellite parallax will be calculated correctly
    for the model evaluated at the data points, but satellite parallax
    is not implemented for the model alone.
    """

    def __init__(self, parameters=None,
                 t_0=None, u_0=None, t_E=None, rho=None, s=None, q=None,
                 alpha=None,
                 pi_E=None, pi_E_N=None, pi_E_E=None,
                 pi_E_ref=None, t_0_par=None, 
                 coords=None, ra=None, dec=None):
        """
        Two ways to define the model:
        1. parameters = a ModelParameters() object
        2. specify t_0, u_0, t_E (optionally: rho, s, q, alpha, pi_E, t_0_par)

        When defining event coordinates, may specify coords as an
        astropy.coordinates.SkyCoord object, otherwise assumes RA is
        in hour angle and DEC is in degrees.
        """
        # Initialize the parameters of the model
        if isinstance(parameters, ModelParameters):
            self._parameters = parameters
        elif parameters is None:
            self._parameters = ModelParameters()
        else:
            raise TypeError(
                "If specified, parameters must be a ModelParameters object.")

        # Set each model parameter
        if t_0 is not None:
            self.t_0 = t_0
        if u_0 is not None:
            self.u_0 = u_0
        if t_E is not None:
            self.t_E = t_E
        if rho is not None:
            self.rho = rho
        self.t_0_par = t_0_par
        
        # Set the parallax
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

        # Set the coordinates of the event
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

        # Set some defaults
        self._parallax = {'earth_orbital':False, 
                          'satellite':False, 
                          'topocentric':False}
        self._satellite_coords = None
        #self._delta_annual = {}
        #self._delta_satellite = {}

    def __repr__(self):
        return '{0}'.format(self._parameters)

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

    @property
    def t_E(self):
        """
        The Einstein timescale. An astropy.Quantity. "day" is the default unit.
        """
        return self._parameters.t_E

    @t_E.setter
    def t_E(self, value):
        self._parameters.t_E = value
        
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

    @property
    def pi_E_N(self):
        """
        The North component of the microlens parallax vector.
        """
        return self._parameters.pi_E_N

    @pi_E_N.setter
    def pi_E_N(self, value):
        self._parameters.pi_E_N = value

    @property
    def pi_E_E(self):
        """
        The East component of the microlens parallax vector.
        """
        return self._parameters.pi_E_E

    @pi_E_E.setter
    def pi_E_E(self, value):
        self._parameters.pi_E_E = value

    @property
    def t_0_par(self):
        """reference time for parameters, in particular microlensing parallax"""
        return self._t_0_par

    @t_0_par.setter
    def t_0_par(self, value):
        self._t_0_par = value

    @property
    def parameters(
        self, t_0=0., u_0=None, t_E=1., rho=None, s=None, q=None, alpha=None, 
        pi_E=None, pi_E_N=None, pi_E_E=None, pi_E_ref=None):
        if u_0 is None:
            return "{0}".format(self._parameters)
        else:
            self._parameters = ModelParameters(
                t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, s=s, q=q, alpha=alpha, 
                pi_E=pi_E, pi_E_N=pi_E_N, pi_E_E=pi_E_E, pi_E_ref=pi_E_ref)

    def magnification(self, time, satellite_coords=None):
        """
        calculate the model magnification for the given time(s).
        """
        if satellite_coords is None:
            satellite_coords = self._satellite_coords

        magnification_curve = MagnificationCurve(
            time, parameters=self._parameters, 
            parallax=self._parallax, t_0_par=self.t_0_par,
            coords=self._coords, 
            satellite_coords=satellite_coords)
        return magnification_curve.magnification

    @property
    def data_magnification(self):
        """a list of magnifications calculated for every dataset time vector"""
        if self._data_magnification is not None:
            return self._data_magnification

        self._data_magnification = []

        for dataset in self._datasets:
            magnification = self.get_data_magnification(dataset)
            self._data_magnification.append(magnification)

        return self._data_magnification
        
    def get_data_magnification(self, dataset):
        """
        Get the model magnification for a given dataset.
        """
        if dataset.is_satellite:
            dataset_satellite_coords = dataset.satellite_skycoord
        else:
            dataset_satellite_coords = None
            
        magnification = self.magnification(
            dataset.time, satellite_coords=dataset_satellite_coords)
        return magnification
        

    def set_datasets(self, datasets):
        """set _datasets property"""
        self._datasets = datasets
        self._data_magnification = None

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

    def parallax(
        self, earth_orbital=False, satellite=False, topocentric=False):
        """
        specifies which types of the parallax will be included in calculations
        """
        self._parallax['earth_orbital'] = earth_orbital
        self._parallax['satellite'] = satellite
        self._parallax['topocentric'] = topocentric


    def plot_magnification(
        self, times=None, t_range=None, t_start=None, t_stop=None, dt=None, 
        n_epochs=None, **kwargs):
        """
        plot the model magnification curve.
        """
        if times is None:
            times = self.set_times(
                parameters=self._parameters, t_range=t_range, t_start=t_start, 
                t_stop=t_stop, dt=dt, 
                n_epochs=n_epochs)

        pl.plot(times, self.magnification(times),**kwargs)
        pl.ylabel('Magnification')
        pl.xlabel('Time')

    def plot_lc(
        self, times=None, t_range=None, t_start=None, t_stop=None, dt=None, 
        n_epochs=None, data_ref=None, f_source=None, f_blend=None, **kwargs):
        """
        plot the model light curve in magnitudes.
        """
        if times is None:
            times = self.set_times(
                parameters=self._parameters, t_range=t_range, t_start=t_start, 
                t_stop=t_stop, dt=dt, 
                n_epochs=n_epochs)

        if (f_source is None) and (f_blend is None):
            f_source, f_blend = self.get_ref_fluxes(data_ref=data_ref)
        elif (f_source is None) or (f_blend is None):
            raise AttributeError(
                'If f_source is set, f_blend must also be set and vice versa.')
            
        flux = f_source * self.magnification(times) + f_blend

        pl.plot(times, Utils.get_mag_from_flux(flux),**kwargs)
        pl.ylabel('Magnitude')
        pl.xlabel('Time')
        
        ymin, ymax = pl.gca().get_ylim()
        if ymax > ymin:
            pl.gca().invert_yaxis()

    def get_ref_fluxes(self, data_ref=None):
        """
        Determine the reference flux system from the
        datasets. data_ref may either be a dataset or the index of a
        dataset (if model.set_datasets() was previously called). If
        data_ref is not set, it will use the first dataset. If you
        call this without calling set_datasets() first, there will be
        an exception and that's on you.
        """
        if data_ref is None:
            data = self._datasets[0]
        elif isinstance(data_ref, MulensData):
            data = data_ref
        else:
            data = self._datasets[data_ref]

        fit = Fit(data=data, magnification=[self.get_data_magnification(data)])
        fit.fit_fluxes()
        f_source = fit._flux_sources[data]
        f_blend = fit._flux_blending[data]

        return f_source, f_blend

    def plot_data(self, data_ref=None):
        """
        Plot the data scaled to the model. If data_ref is not
        specified, uses the first dataset as the flux
        reference. 
        """

        all_data_fit = Fit(
            data=self._datasets, magnification=self.data_magnification())

    def set_times(
        self, parameters=None, t_range=None, t_start=None, t_stop=None, 
        dt=None, n_epochs=None):
        """
        If given, set up a time vector based on t_start, t_stop,
        and (dt or n_epochs). If not given, intialize the time
        vector based on the model parameters.
        """
            #initialize t_start, t_stop, dt if not set
        if t_range is not None:
            t_start = t_range[0]
            t_stop = t_range[1]

        n_tE = 1.5
        if t_start is None:
            t_start = self.t_0 - (n_tE * self.t_E)
        if t_stop is None:
            t_stop = self.t_0 + (n_tE * self.t_E)
                
        if dt is None:
            if n_epochs is None:
                n_epochs = 1000
            dt = (t_stop - t_start) / float(n_epochs)

        return np.arange(t_start, t_stop+dt, dt)
