import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import GeocentricTrueEcliptic
import matplotlib.pyplot as pl

from MulensModel.modelparameters import ModelParameters
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.trajectory import Trajectory
from MulensModel.caustics import Caustics
from MulensModel.mulensparallaxvector import MulensParallaxVector
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs

#JCY: some probable problems with annual parallax due to interface
#with astropy functions get_body_barycentric and get_jd12. These must
#depend on both reference frame and time standard. But it is only
#possible to set time standard and jd vs. mjd.



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

    def __init__(
        self, parameters=None, t_0=None, u_0=None, t_E=None, rho=None, s=None, 
        q=None, alpha=None, pi_E=None, pi_E_N=None, pi_E_E=None, pi_E_ref=None, 
        t_0_par=None, coords=None, ra=None, dec=None, ephemerides_file=None):
        """
        Two ways to define the model:
        1. parameters = a ModelParameters() object
        2. specify t_0, u_0, t_E (optionally: rho, s, q, alpha, pi_E, t_0_par)

        When defining event coordinates, may specify coords as an
        astropy.coordinates.SkyCoord object, otherwise assumes RA is
        in hour angle and DEC is in degrees.

        Default values for parallax are all True. Use model.parallax() to turn
        different parallax effects ON/OFF.
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
        if s is not None:
            self.s = s

        if q is not None:
            if isinstance(q, (list, np.ndarray)):
                if len(q) > 1:
                    raise NotImplementedError(
                        'Too many q. Does not support more than 2 bodies.')
                else:
                    q = q[0]
            self.q = q

        if alpha is not None:
            self.alpha = alpha
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

        self.ephemerides_file = ephemerides_file
        self._satellite_skycoord = None
        
        # Set some defaults
        self._parallax = {'earth_orbital':True, 
                          'satellite':True, 
                          'topocentric':True}
        self._default_magnification_method = 'point_source'
        self._methods = None
        self.caustics = None

        #Set dictionary to store plotting properties
        self.reset_plot_properties()
        
        self._limb_darkening_coeffs = LimbDarkeningCoeffs()

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
    def s(self):
        """lens components separation in units of theta_E"""
        return self._parameters.s

    @s.setter
    def s(self, value):
        self._parameters.s = value

    @property
    def q(self):
        """mass ratio of lens components"""
        return self._parameters.q

    @q.setter
    def q(self, value):
        self._parameters.q = value

    @property
    def alpha(self):
        """angle between lens axis and source trajectory"""
        return self._parameters.alpha

    @alpha.setter
    def alpha(self, value):
        self._parameters.alpha = value
    
    @property
    def parameters(self):
        """
        The parameters of the model. A Model Parameters object.
        """
        return self._parameters

    def set_parameters(
        self, t_0=0., u_0=None, t_E=1., rho=None, s=None, q=None, alpha=None, 
        pi_E=None, pi_E_N=None, pi_E_E=None, pi_E_ref=None):
        """
        Set the parameters of the model. Any parameter not explicitly
        specified will be set to None.
        """
        self._parameters = ModelParameters(
            t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, s=s, q=q, alpha=alpha, 
            pi_E=pi_E, pi_E_N=pi_E_N, pi_E_E=pi_E_E, pi_E_ref=pi_E_ref)
            
    def get_satellite_coords(self,times):
        """
        get satellite SkyCoords for the given times
        """
        if self.ephemerides_file is None:
            return None
        else:
            satellite_skycoords = SatelliteSkyCoord(
                 ephemerides_file=self.ephemerides_file)
            return satellite_skycoords.get_satellite_coords(times)

    def magnification(self, time, satellite_skycoord=None, gamma=0.):
        """
        calculate the model magnification for the given time(s).
        """
        #Check for type
        if not isinstance(time, np.ndarray):
            if isinstance(time, (np.float, float)):
                time = np.array([time])
            elif isinstance(time, list):
                time = np.array(time)
            else:
                raise TypeError('time must be a float, list, or np.ndarray')

        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(time)

        magnification_curve = MagnificationCurve(
            time, parameters=self._parameters, 
            parallax=self._parallax, t_0_par=self.t_0_par,
            coords=self._coords, 
            satellite_skycoord=satellite_skycoord,
            gamma=gamma)
        magnification_curve.set_magnification_methods(self._methods, 
                                        self._default_magnification_method)
        return magnification_curve.magnification

    @property
    def data_magnification(self):
        """a list of magnifications calculated for every dataset time vector"""
        self._data_magnification = []

        for dataset in self._datasets:
            magnification = self.get_data_magnification(dataset)
            self._data_magnification.append(magnification)

        return self._data_magnification
        
    def get_data_magnification(self, dataset):
        """
        Get the model magnification for a given dataset.
        """
        if dataset.ephemerides_file is not None:
            dataset_satellite_skycoord = dataset.satellite_skycoord
        else:
            dataset_satellite_skycoord = None
            
        if dataset.bandpass is None:
            gamma = 0.
        else:
            gamma = self._limb_darkening_coeffs.limb_coef_gamma(dataset.bandpass)
            
        magnification = self.magnification(
            dataset.time, satellite_skycoord=dataset_satellite_skycoord, gamma=gamma)
        return magnification
        
    def set_datasets(self, datasets, data_ref=0):
        """set _datasets property"""
        self._datasets = datasets
        self._data_magnification = None
        self.data_ref = data_ref

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
        if l > 180. * u.deg:
            l = l - 360. * u.deg
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
        self, earth_orbital=None, satellite=None, topocentric=None):
        """
        specifies which types of the parallax will be included in
        calculations. Three kinds of effects are allowed:
        earth_orbital - the motion of the Earth about the Sun
        satellite - difference due to the separation between the Earth
            and a satellite (changes as a function of time)
        topocentric - difference due to the separation between two
            observatories on the Earth.
        """
        if earth_orbital is None and satellite is None and topocentric is None:
            return self._parallax
        else:
            if earth_orbital is not None:
                self._parallax['earth_orbital'] = earth_orbital
            if satellite is not None:
                self._parallax['satellite'] = satellite
            if topocentric is not None:
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

        pl.plot(times, self.magnification(times), **kwargs)
        pl.ylabel('Magnification')
        pl.xlabel('Time')

    def plot_lc(
        self, times=None, t_range=None, t_start=None, t_stop=None, dt=None, 
        n_epochs=None, data_ref=None, f_source=None, f_blend=None, 
        **kwargs):
        """
        plot the model light curve in magnitudes. See get_ref_fluxes
        for details of data_ref.
        """
        if times is None:
            times = self.set_times(
                parameters=self._parameters, t_range=t_range, t_start=t_start, 
                t_stop=t_stop, dt=dt, 
                n_epochs=n_epochs)

        if data_ref is not None:
            self.data_ref = data_ref

        if (f_source is None) and (f_blend is None):
            (f_source, f_blend) = self.get_ref_fluxes(data_ref=data_ref)
        elif (f_source is None) or (f_blend is None):
            raise AttributeError(
                'If f_source is set, f_blend must also be set and vice versa.')
            
        flux = f_source * self.magnification(times) + f_blend

        pl.plot(times, Utils.get_mag_from_flux(flux), **kwargs)
        pl.ylabel('Magnitude')
        pl.xlabel('Time')
        
        (ymin, ymax) = pl.gca().get_ylim()
        if ymax > ymin:
            pl.gca().invert_yaxis()

    def get_ref_fluxes(self, data_ref=None):
        """
        Determine the reference flux system from the
        datasets. data_ref may either be a dataset or the index of a
        dataset (if Model.set_datasets() was previously called). If
        data_ref is not set, it will use the first dataset. If you
        call this without calling set_datasets() first, there will be
        an exception and that's on you.
        """
        if data_ref is None:
            if isinstance(self.data_ref, MulensData):
                data = self.data_ref
            else:
                data = self._datasets[self.data_ref]
        elif isinstance(data_ref, MulensData):
            data = data_ref
            self.data_ref = data_ref
        else:
            data = self._datasets[data_ref]
            self.data_ref = data_ref

        fit = Fit(data=data, magnification=[self.get_data_magnification(data)])
        fit.fit_fluxes()
        f_source = fit.flux_of_sources(data)
        f_blend = fit.blending_flux(data)

        return (f_source, f_blend)
        
    def reset_plot_properties(self):
        """resets internal settings for plotting"""
        self.plot_properties = {}

    def _store_plot_properties(
        self, color_list=None, marker_list=None, size_list=None, 
        label_list=None, **kwargs):
        """
        Store plot properties for each data set.
        """
        if color_list is not None:
            self.plot_properties['color_list'] = color_list
        if marker_list is not None:
            self.plot_properties['marker_list'] = marker_list
        if size_list is not None:
            self.plot_properties['size_list'] = size_list
        if label_list is not None:
            self.plot_properties['label_list'] = label_list
        if len(kwargs) > 0:
            self.plot_properties['other_kwargs'] = kwargs
    
    def _set_plot_kwargs(self, index, show_errorbars=True):
        """
        Set kwargs arguments for plotting. If set, use previous values. But 
        new values take precedence. 
        
        Automatically handles (some) differences in keywords for pl.errorbar 
        vs. pl.scatter: fmt/marker, markersize/s
        """                
        #Set different keywords for pl.errorbar vs. pl.scatter
        if show_errorbars:
            marker_key = 'fmt'
            size_key = 'markersize'
        else:
            marker_key = 'marker'
            size_key = 's'
    
        #Create new kwargs dictionary
        new_kwargs = {}  
        
        #Set defaults
        if index == 0:
            pl.gca().set_color_cycle(None)
        new_kwargs[marker_key] = 'o'
        new_kwargs[size_key] = 3
        
        #Set custom
        if len(self.plot_properties) > 0:
            if 'color_list' in self.plot_properties.keys():
                new_kwargs['color'] = self.plot_properties['color_list'][index]

            if 'marker_list' in self.plot_properties.keys():
                new_kwargs[marker_key] = self.plot_properties['marker_list'][index]

            if 'size_list' in self.plot_properties.keys():
                new_kwargs[size_key] = self.plot_properties['size_list'][index]

            if 'label_list' in self.plot_properties.keys():
                new_kwargs['label'] = self.plot_properties['label_list'][index]
                
            if 'other_kwargs' in self.plot_properties.keys():
                for (key, value) in self.plot_properties['other_kwargs'].items():
                    if key == 'markersize' or key == 's':
                        new_kwargs[size_key] = value
                    elif key == 'marker' or key == 'fmt':
                        new_kwargs[marker_key] = value
                    else:
                        new_kwargs[key] = value
                        
        return new_kwargs

    def plot_data(
        self, data_ref=None, show_errorbars=True, color_list=None,
        marker_list=None, size_list=None, label_list=None, **kwargs):
        """
        Plot the data scaled to the model. If data_ref is not
        specified, uses the first dataset as the reference for flux scale. 

        If show_errorbars is True (default), plots with matplotlib.errorbar(). 
        If show_errorbars is False, plots with matplotib.scatter(). 
        
        Allows for different point types for each dataset. These may be set
        using color_list, marker_list, and size_list. May also use **kwargs
        or some combination of the lists and **kwargs. e.g. set color_list to 
        specify which color each data set should be plotted in, but use 
        fmt='s' to make all data points plotted as squares.
        
        Automatically handles some keyword variations in errorbar() vs. 
        scatter(): e.g. fmt/marker, markersize/s (see _set_plot_kwargs),
        
        **kwargs (and point type lists) are remembered and used in subsequent 
        calls to both plot_data() and plot_residuals(). 
        """
        if data_ref is not None:
            self.data_ref = data_ref

        self._store_plot_properties(
            color_list=color_list, marker_list=marker_list, 
            size_list=size_list, label_list=label_list,
            **kwargs)

        #Reference flux scale
        (f_source_0, f_blend_0) = self.get_ref_fluxes(data_ref=data_ref)

        #Get fluxes for all datasets
        fit = Fit(
            data=self._datasets, magnification=self.data_magnification)
        fit.fit_fluxes()
        
        #plot defaults
        t_min = 3000000.
        t_max = 0.
        
        #plot each dataset
        for (i, data) in enumerate(self._datasets):
            #Calculate scaled flux
            f_source = fit.flux_of_sources(data)
            f_blend = fit.blending_flux(data)
            flux = f_source_0 * (data.flux - f_blend) / f_source + f_blend_0

            new_kwargs = self._set_plot_kwargs(i, show_errorbars=show_errorbars)
            #Plot
            if show_errorbars:
                err_flux = f_source_0 * data.err_flux / f_source
                (mag, err) = Utils.get_mag_and_err_from_flux(flux, err_flux)
                pl.errorbar(data.time, mag, yerr=err, **new_kwargs) 
                
            else:
                mag = Utils.get_mag_from_flux(flux)
                pl.scatter(data.time, mag, lw=0., **new_kwargs)

            #Set plot limits
            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))

        #Plot properties
        pl.ylabel('Magnitude')
        pl.xlabel('Time')
        pl.xlim(t_min, t_max)

        (ymin, ymax) = pl.gca().get_ylim()
        if ymax > ymin:
            pl.gca().invert_yaxis()

    def plot_residuals(
        self, show_errorbars=True, color_list=None, marker_list=None, 
        size_list=None, label_list=None, data_ref=None, **kwargs):
        """
        Plot the residuals (in magnitudes) of the model. Uses the best f_source,
        f_blend for each dataset (not scaled to a particular 
        photometric system).

        For explanation of **kwargs, and also [var]_list see doctrings in 
        plot_data(). 
        """
        if data_ref is not None:
            self.data_ref = data_ref

        self._store_plot_properties(
            color_list=color_list, marker_list=marker_list, 
            size_list=size_list, label_list=label_list,
            **kwargs)
            
        #Get fluxes for all datasets
        fit = Fit(
            data=self._datasets, magnification=self.data_magnification)
        fit.fit_fluxes()

        #Plot limit parameters
        delta_mag = 0.
        t_min = 3000000.
        t_max = 0.

        #Plot zeropoint line
        pl.plot([0., 3000000.], [0., 0.], color='black')
        
        #Plot residuals
        for (i, data) in enumerate(self._datasets):
            #Calculate model magnitude
            f_source = fit.flux_of_sources(data)
            f_blend = fit.blending_flux(data)
            model_flux = f_source * self.magnification(data.time) + f_blend
            model_mag = Utils.get_mag_from_flux(model_flux)

            #Calculate Residuals
            (mag, err) = Utils.get_mag_and_err_from_flux(data.flux, 
                                                       data.err_flux)
            residuals = model_mag - mag
            delta_mag = max(delta_mag, np.max(np.abs(residuals)))

            #Plot
            new_kwargs = self._set_plot_kwargs(i, show_errorbars=show_errorbars)
            if show_errorbars:
                pl.errorbar(data.time, residuals, yerr=err, 
                            **new_kwargs) 
            else:
                pl.scatter(data.time, residuals, lw=0, **new_kwargs)

            #Set plot limits
            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))
       
        if delta_mag > 1.:
            delta_mag = 0.5

        #Plot properties
        pl.ylim(-delta_mag, delta_mag)
        pl.xlim(t_min, t_max)
        pl.ylabel('Residuals')
        pl.xlabel('Time')

    def plot_trajectory(
        self, times=None, t_range=None, t_start=None, t_stop=None, dt=None, 
        n_epochs=None, caustics=False, show_data=False, arrow=True,
        satellite_skycoord=None,**kwargs):
        """
        Plot the source trajectory.

        Optional keyword arguments:

          times, t_range, t_start, t_stop, dt, n_epochs may all be
          used to specify exactly when to plot the source
          trajectory. times=specific dates, t_range=range of times,
          (t_start, t_stop)=range of times with optional dt OR
          n_epochs.

          caustics = plot the caustic structure in addition to the
          source trajectory. default=False (off). For finer control of
          plotting features, e.g. color, use self.plot_caustics()
          instead.

          show_data = mark epochs of data (Not Implemented, marker
          types should match data plotting.)

          arrow = show the direction of the motion. default=True (on)

          satellite_skycoord should allow user to specify the trajectory
          is calculated for a satellite. (Not checked)

          **kwargs controls plotting features of the trajectory.
        """
        if times is None:
            times = self.set_times(
                parameters=self._parameters, t_range=t_range, t_start=t_start, 
                t_stop=t_stop, dt=dt, 
                n_epochs=n_epochs)

        if satellite_skycoord is None:
            satellite_skycoord = self.get_satellite_coords(times)

        trajectory = Trajectory(
            times, parameters=self._parameters, parallax=self._parallax, 
            t_0_par=self.t_0_par, coords=self._coords, 
            satellite_skycoord=satellite_skycoord)

        pl.plot(trajectory.x, trajectory.y, **kwargs)
        
        if arrow:
            index = len(times)/2
            pl.scatter(
                trajectory.x[index], trajectory.y[index], 
                marker=(3, 0, self.alpha), s=50)

        if caustics:
            self.plot_caustics(marker='.', color='red')

    def plot_caustics(self, n_points=5000, **kwargs):
        """
        Plot the caustic structure. n_points specifies the number of
        points to generate in the caustic.
        """
        if self.caustics is None:
            self.caustics = Caustics(q=self.q, s=self.s)

        self.caustics.plot(n_points=n_points,**kwargs)
        
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

    def set_default_magnification_method(self, method):
        """stores information on method to be used, when no metod is
        directly specified"""
        self._default_magnification_method = method

    def set_magnification_methods(self, methods):
        """sets methods used for magnification calculation
        
        Parameter method is a list that contains epochs and names of methods
        to be used:
        methods = [2455746., 'Quadrupole', 2455746.6, 'Hexadecapole', 
                   2455746.7, 'VBBL', 2455747., 'Hexadecapole', 2455747.15, 
                   'Quadrupole', 2455748.]
        """
        self._methods = methods

    def set_limb_coef_gamma(self, bandpass, coef):
        """store gamma LD coef for given band"""
        self._limb_darkening_coeffs.set_limb_coef_gamma(bandpass, coef)

    def set_limb_coef_u(self, bandpass, coef):
        """store u LD coef for given band"""
        self._limb_darkening_coeffs.set_limb_coef_u(bandpass, coef)

    def limb_coef_gamma(self, bandpass):
        """get gamma LD coef for given band"""
        self._limb_darkening_coeffs.limb_coef_gamma(bandpass)

    def limb_coef_u(self, bandpass):
        """get u LD coef for given band"""
        self._limb_darkening_coeffs.limb_coef_u(bandpass)

    @property
    def bandpasses(self):
        """ """
        pass
