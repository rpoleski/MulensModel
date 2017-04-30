import numpy as np
from scipy.special import ellipeinc # This is incomplete elliptic integral of the second kind.

from MulensModel.trajectory import Trajectory
from MulensModel.binarylens import BinaryLens
from MulensModel.modelparameters import ModelParameters


class MagnificationCurve(object):
    """
    The magnification curve calculated from the model light curve.
    """
    def __init__(
        self, times, parameters=None, parallax=None, t_0_par=None,
        coords=None, satellite_skycoord=None):
        """
        Required arguments: 
           times - the times at which to generate the magnification curve, e.g. a vector.
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

        #Calculate the source trajectory (i.e. u(t))
        self.trajectory = Trajectory(
            self.times, parameters=parameters, parallax=parallax, 
            t_0_par=t_0_par, coords=coords, satellite_skycoord=satellite_skycoord)

        #Initialize the magnification vector
        self._magnification = None
            
    @property
    def magnification(self):
        """provide vector of magnifications"""
        return self.get_magnification()

    def get_magnification(self):
        """calculate magnification"""
        if self.parameters.n_lenses == 1:
            magnification = self.get_point_lens_magnification()
        elif self.parameters.n_lenses == 2:
            magnification = self.get_binary_lens_magnification()
        else:
            raise Exception(
                "magnification for more than 2 lenses not handled yet")
        self._magnification = magnification
        return self._magnification

    def get_point_lens_magnification(self):
        """Calculate the Point Lens magnification. """
        u2 = (self.trajectory.x**2 + self.trajectory.y**2)
        return (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    
    def _get_point_lens_finite_source_magnification(self, rho, mask):
        """calculate magnification for point lens and finite source. 
        Variable mask defines which epochs to use
        
        The approximation was propsed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G
        
        and later the integral calculation was simplified by:
        
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y
        """

        u2 = (self.trajectory.x[mask]**2 + self.trajectory.y[mask]**2)
        pspl_magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
        z = np.sqrt(u2) / rho
        
        k = np.ones_like(z)
        for (i, value) in enumerate(z):
            if value > 1.:
                k[i] = 1. / value
                
        B_0 = 4. * z * ellipeinc(k, z) / np.pi # I'm not sure if the order of arguments is correct. This can be easily checked.
        
        magnification = pspl_magnification * B_0
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
        return magnification

    def get_binary_lens_magnification(self):
        """Calculate the Binary magnification. """
        #Set up the Binary lens system
        q = self.parameters.q
        m_1 = 1. / (1. + q)
        m_2 = q / (1. + q)
        binary_lens = BinaryLens(mass_1=m_1, mass_2=m_2, 
                                    separation=self.parameters.s)
        
        #Calculate the magnification
        magnification = []
        for i in range(len(self.trajectory.x)):
            x = self.trajectory.x[i]
            y = self.trajectory.y[i]
            m = binary_lens.point_source_magnification(source_x=x, source_y=y)
            magnification.append(m)
            
        return np.array(magnification)
        