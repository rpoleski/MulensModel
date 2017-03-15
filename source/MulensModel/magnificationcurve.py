import numpy as np

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
    
    def get_binary_lens_magnification(self):
        """Calculate the Binary magnification. """
        #Set up the Binary lens system
        q = self.parameters.q
        m1 = 1. / (1. + q)
        m2 = q / (1. + q)
        binary_lens = BinaryLens(mass_1=m1, mass_2=m2, 
                                    separation=self.parameters.s)
        
        #Calculate the magnification
        magnification = []
        for i in range(len(self.trajectory.x)):
            x = self.trajectory.x[i]
            y = self.trajectory.y[i]
            m = binary_lens.point_source_magnification(source_x=x, source_y=y)
            magnification.append(m)
            
        return np.array(magnification)
        
