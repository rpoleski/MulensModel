import numpy as np
from MulensModel.trajectory import Trajectory
from MulensModel.binarylensequation import BinaryLensEquation

class MagnificationCurve(object):
    """
    The magnification curve calculated from the model light curve.
    """
    def __init__(
        self, times, parameters=None, parallax=None, t_0_par=None,
        coords=None, satellite_coords=None):

        self.times = times
        self.parameters = parameters
        self.trajectory = Trajectory(
            self.times, parameters=parameters, parallax=parallax, 
            t_0_par=t_0_par, coords=coords, satellite_coords=satellite_coords)
        self.get_magnification()


    def get_magnification(self):
        if self.parameters.n_lenses == 1:
            magnification = self.get_point_lens_magnification()
        elif self.parameters.n_lenses == 2:
            magnification = self.get_binary_lens_magnification()
        else:
            raise Exception(
                "magnification for more than 2 lenses not handled yet")
        self.magnification = magnification
        return self.magnification


    def get_point_lens_magnification(self):
        """Calculate the Point Lens magnification. """
        u2 = (self.trajectory.x**2 + self.trajectory.y**2)
        return (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    

    def get_binary_lens_magnification(self):
        """Calculate the Binary magnification. """
        q = self.parameters.q
        m1 = 1. / (1. + q)
        m2 = q / (1. + q)
        binary_lens_eq = BinaryLensEquation(
            mass_1=m1, mass_2=m2, separation=self.parameters.s, 
            source_x=self.trajectory.x, 
            source_y=self.trajectory.y)
        
        return binary_lens_eq.total_magnification
