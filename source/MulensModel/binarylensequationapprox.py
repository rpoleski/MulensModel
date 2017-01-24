import numpy as np

from MulensModel.BinaryLensEquation import BinaryLensEquation


class BinaryLensApprox(object):
    """hexadecpole and quadrupole approximation of binary-lens/finite source 
    calculations - based on Gould 2008 ApJ 681, 1593"""
    def __init__(self, mass_1=None, mass_2=None, separation=None, rho=None, 
                 source_x=None, source_y=None, gamma=None, 
                 degree='hexadecpole'): 
        self.mass_1 = mass_1
        self.mass_2 = mass_2
        self.separation = separation
        self.rho = rho
        self.source_x = source_x
        self.source_y = source_y
        if degree in ['hexadecpole', 'quadrupole']:
            self.degree = degree
        else:
            raise ValueError('Unknown degree in binary lens approximation')
        
    def _get_a_w_plus(self, radius):
        """This is Gould 2008 eq. 7"""
        x = self.source_x + radius
        y = self.source_y
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_1 = binary_lens_eq.total_magnification
        x = self.source_x
        y = self.source_y + radius
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_2 = binary_lens_eq.total_magnification
        x = self.source_x - radius
        y = self.source_y 
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_3 = binary_lens_eq.total_magnification
        x = self.source_x
        y = self.source_y - radius
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_4 = binary_lens_eq.total_magnification
        # Yes, above we have almost the same block of code repeated 4 times. 
        return 0.25 * (a_1 + a_2 + a_3 + a_4) - self.value_center

        
    def _get_a_w_times(self, radius):
        """This is Gould 2008 eq. 8"""
        shift = radius / 2**0.5
        x = self.source_x + shift
        y = self.source_y + shift
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_1 = binary_lens_eq.total_magnification
        x = self.source_x - shift
        y = self.source_y + shift
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_2 = binary_lens_eq.total_magnification
        x = self.source_x - shift
        y = self.source_y - shift 
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_3 = binary_lens_eq.total_magnification
        x = self.source_x + shift
        y = self.source_y - shift
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=x, source_y=y)
        a_4 = binary_lens_eq.total_magnification
        # Yes, above we have almost the same block of code repeated 4 times. 
        return 0.25 * (a_1 + a_2 + a_3 + a_4) - self.value_center


    def get_magnification(self):
        """calculates approximated finite source magnification"""
        binary_lens_eq = BinaryLensEquation(mass_1=self.mass_1, 
                mass_2=self.mass_2, separation=self.separation, 
                source_x=self.source_x, source_y=self.source_y)
        self.value_center = binary_lens_eq.total_magnification
        a_rho_half_plus = self._get_a_w_plus(radius=0.5*self.rho)
        a_rho_plus = self._get_a_w_plus(radius=self.rho)
        
        # This is Gould 2008 eq. 9:
        a_2_rho_sq = (16. * a_rho_half_plus - a_rho_plus) / 3. 
        
        # Gould 2008 eq. 6 (part 1/2):
        a_finite = self.value_center + a_2_rho_sq * (1. - 0.2 * self.gamma) 
        if self.degree == 'quadrupole':
            return a_finite
            
        a_rho_times = self._get_a_w_times(radius=self.rho)
        # This is Gould 2008 eq. 9:
        a_4_rho_pow4 = 0.5 * (a_rho_plus + a_rho_times) - a_2_rho_sq
        a_finite += a_4_rho_pow4 * (1. - 11. * self.gamma / 35.)
        return a_finite
        