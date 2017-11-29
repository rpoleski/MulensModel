import numpy as np
import warnings
from scipy.special import ellipe # This is incomplete elliptic integral of the second kind.
from scipy import integrate

from MulensModel.trajectory import Trajectory
from MulensModel.binarylens import BinaryLens
from MulensModel.modelparameters import ModelParameters


class MagnificationCurve(object):
    """
    The magnification curve calculated from the model light curve.
    """
    def __init__(self, times, parameters=None, parallax=None,  
                    coords=None, satellite_skycoord=None, gamma=0.):
        """
        Required arguments: 
           times - the times at which to generate the magnification curve, e.g. a vector.
           parameters - a ModelParameters object specifying the microlensing parameters

        Optional parallax keywords:
           parallax - boolean dictionary specifying what parallax effects should be used.
           coords - sky coordinates of the event
           satellite_skycoord - sky coordinates of the satellite
               specified by the ephemrides file. see
               MulensData.satellite_skycoord.
           gamma - limb darkening coefficient in gamma convention
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

        #Calculate the source trajectory (i.e. u(t))
        self.trajectory = Trajectory(
            self.times, parameters=parameters, parallax=parallax, 
            coords=coords, satellite_skycoord=satellite_skycoord)

        #Initialize the magnification vector
        self._magnification = None

        #Set methods' variables:
        self._methods_epochs = None
        self._methods_names = None
        self._default_magnification_method = None

        self._gamma = gamma

    def set_magnification_methods(self, methods, default_method):
        """sets methods used for magnification calculation;
        epochs is a numpy ndarray of n epochs that specify when (n-1) 
        methods will be used"""
        self._default_method = default_method
        if methods is None:
            self._methods_epochs = None
            self._methods_names = None
            return
        
        if not isinstance(methods, list):
            msg = 'MagnificationCurve.set_magnification_methods() requires list as a parameter'
            raise TypeError(msg)
        epochs = methods[0::2]
        names = methods[1::2]
        
        for epoch in epochs:
            if not isinstance(epoch, float):
                raise TypeError('Wrong epoch: {:}'.format(epoch))
        for method in names:
            if not isinstance(method, str):
                raise TypeError('Wrong method: {:}'.format(method))
        for (e_beg, e_end) in zip(epochs[::2], epochs[1::2]):
            if e_beg >= e_end:
                msg = 'Incorrect epochs provided: {:} and {:} (first should be earlier)'
                raise ValueError(msg.format(e_beg, e_end))
        
        self._methods_epochs = np.array(epochs)
        self._methods_names = names        
        self._default_method = default_method

    @property
    def magnification(self):
        """provide vector of magnifications"""
        return self.get_magnification()
        # THIS HAS TO BE REWRITTEN - USE LAZY LOADING! (here or in model.py)

    def get_magnification(self):
        """Calculate magnification.
        
        Returns :
            magnification: *np.ndarray*
                Vector of magnifications. 
        """
        if self.parameters.n_lenses == 1:
            magnification = self.get_point_lens_magnification()
        elif self.parameters.n_lenses == 2:
            magnification = self.get_binary_lens_magnification()
        else:
            raise NotImplementedError(
                "magnification for more than 2 lenses not handled yet")
        self._magnification = magnification
        return self._magnification

    def get_point_lens_magnification(self):
        """Calculate the Point Lens magnification. 
        
        Returns :
            magnification: *np.ndarray*
                Vector of magnifications.
        """

        if self.parameters.rho is not None:
            if self._methods_epochs is None:
                warnings.warn('1 rho set but no finite-source method is set')
            elif set(self._methods_for_epochs()) == set(['point_source']):
                warnings.warn('2 rho set but no finite-source method is set')

        u2 = (self.trajectory.x**2 + self.trajectory.y**2)
        # This is Paczynski equation, i.e., point-source/point-lens (PSPL) 
        # magnification:
        pspl_magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
        if self._methods_epochs is None:
            return pspl_magnification
            
        magnification = pspl_magnification
        u_all = np.sqrt(u2)
        methods = np.array(self._methods_for_epochs())

        for method in set(methods):
            if method.lower() == 'point_source':
                pass # This cases are already taken care of. 
            elif method.lower() == 'finite_source_uniform_Gould94'.lower():
                selection = (methods == method)
                magnification[selection] = self._get_point_lens_finite_source_magnification(
                    rho=self.parameters.rho, u=u_all[selection], 
                    pspl_magnification=pspl_magnification[selection])
            elif method.lower() == 'finite_source_LD_Gould94'.lower():
                selection = (methods == method)
                magnification[selection] = self._get_point_lens_limb_darkening_magnification(
                    rho=self.parameters.rho, u=u_all[selection], 
                    pspl_magnification=pspl_magnification[selection])
                
            else:
                msg = 'Unknown method specified for single lens: {:}'
                raise ValueError(msg.format(method))        
        
        return magnification       
    
    def _B_0_function(self, z):
        """calculate B_0(z) function defined in:
            
        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G
        
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y
        """
        
        out = 4. * z / np.pi
        function = lambda x: (1.-value**2*np.sin(x)**2)**.5

        for (i, value) in enumerate(z):
            if value < 1.:
                out[i] *= ellipe(value*value)
            else:
                out[i] *= integrate.quad(function, 0., np.arcsin(1./value))[0]
        return out
        
    def _B_1_function(self, z, B_0=None):
        """calculate B_1(z) function defined in:
            
        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G
        
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y"""
        if B_0 is None:
            B_0 = self._B_0_function(z)
        
        function = (lambda r, theta: r * np.sqrt(1.-r**2) / 
            self.parameters.rho / np.sqrt(r**2+zz**2-2.*r*zz*np.cos(theta)))
        lim_0 = lambda x: 0
        lim_1 = lambda x: 1
        W_1 = 0. * z
        for (i, zz) in enumerate(z):
            W_1[i] = integrate.dblquad(function, 0., 2.*np.pi, lim_0, lim_1)[0]
        W_1 /= np.pi
        return B_0 - 1.5 * z * self.parameters.rho * W_1
        
    def _get_point_lens_finite_source_magnification(
                self, rho, u, pspl_magnification):
        """calculate magnification for point lens and finite source.         
        The approximation was propsed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G
        
        and later the integral calculation was simplified by:
        
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y
        """

        z = u / rho
        
        magnification = pspl_magnification * self._B_0_function(z)
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
        return magnification

    def _get_point_lens_limb_darkening_magnification(
                self, rho, u, pspl_magnification):        
        """calculate magnification for point lens and finite source with 
        limb darkening. The approximation was propsed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G
        
        and later the integral calculation was simplified by:
        
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y
        """
        z = u / rho
        B_0 = self._B_0_function(z)
        B_1 = self._B_1_function(z, B_0=B_0)
        magnification = pspl_magnification * (B_0 - self._gamma * B_1)
        return magnification
        
    def get_binary_lens_magnification(self):
        """Calculate the binary lens magnification. 
        
        Returns :
            magnification: *np.ndarray*
                Vector of magnifications.
        """
        #Set up the binary lens system
        q = self.parameters.q
        m_1 = 1. / (1. + q)
        m_2 = q / (1. + q)
        binary_lens = BinaryLens(mass_1=m_1, mass_2=m_2, 
                                    separation=self.parameters.s)
        methods = self._methods_for_epochs()
        
        if self.parameters.rho is not None:
            if self._methods_epochs is None:
                warnings.warn('rho set but no finite-source method is set')
            elif set(methods) != set(['point_source']):
                warnings.warn('rho set but no finite-source method is set')
       
        #Calculate the magnification
        magnification = []        
        for index in range(len(self.times)):
            x = self.trajectory.x[index]
            y = self.trajectory.y[index]
            method = methods[index].lower()
            
            if method == 'point_source':
                m = binary_lens.point_source_magnification(x, y)
            elif method == 'quadrupole':
                m = binary_lens.hexadecapole_magnification(x, y, 
                        rho=self.parameters.rho, quadrupole=True,
                        gamma=self._gamma) 
            elif method == 'hexadecapole':
                m = binary_lens.hexadecapole_magnification(x, y, 
                        rho=self.parameters.rho, 
                        gamma=self._gamma) 
            elif method == 'vbbl':
                m = binary_lens.vbbl_magnification(x, y, 
                        rho=self.parameters.rho, 
                        gamma=self._gamma)
            else:
                msg = 'Unknown method specified for binary lens: {:}'
                raise ValueError(msg.format(method))
            
            magnification.append(m)
            
        return np.array(magnification)

    def _methods_for_epochs(self):
        """for given epochs, decide which methods should be used to calculate
        magnification, but don't run the calculations"""
        out = [self._default_method] * len(self.times)
        if self._methods_epochs is None:
            return out

        brackets = np.searchsorted(self._methods_epochs, self.times)
        n_max = len(self._methods_epochs)
        
        out = [self._methods_names[value-1] 
                    if (value>0 and value<n_max) 
                    else self._default_method 
                    for value in brackets]
        
        return out

