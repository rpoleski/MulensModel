#! /usr/bin/env python

import numpy as np
from math import fsum

from MulensModel.utils import Utils


class BinaryLens(object):
    """binary lens equation - its solutions, images, parities etc."""
    def __init__(self, mass_1=None, mass_2=None, separation=None):
        """mass_1, mass_2, and separation are relative to 
        some mass (and corresponding Einstein radius)"""
        self.mass_1 = mass_1
        self.mass_2 = mass_2
        self.separation = separation
        self._total_mass = None
        self._mass_difference = None
        self._position_z1_WM95 = None
        self._position_z2_WM95 = None
        self._last_polynomial_input = None

    def _calculate_variables(self, source_x, source_y):
        """calculates values of variables needed for polynomial coefficients"""
        self._total_mass = 0.5 * (self.mass_1 + self.mass_2)
        self._mass_difference = 0.5 * (self.mass_2 - self.mass_1)
        self._position_z1_WM95 = -0.5 * self.separation + 0.j
        self._position_z2_WM95 = 0.5 * self.separation + 0.j
        self._zeta_WM95 = source_x - self.separation + source_y * 1.j

    def _get_polynomial_WM95(self, source_x, source_y):
        """calculate coefficients of the polynomial"""
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m
        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff
        pos_z1 = self._position_z1_WM95
        pos_z2 = self._position_z2_WM95
        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2
        zeta = self._zeta_WM95
        zeta_conj = np.conjugate(zeta)
        zeta_conj_pow2 = zeta_conj * zeta_conj

        coef_5 = sum([z1_pow2, -zeta_conj_pow2])
        coef_4 = sum([-2. * total_m * zeta_conj, zeta * zeta_conj_pow2, 
                      -2. * m_diff * pos_z1, -zeta * z1_pow2])
        coef_3 = sum([4. * total_m * zeta * zeta_conj, 
                      4. * m_diff * zeta_conj * pos_z1, 
                      2. * zeta_conj_pow2 * z1_pow2, -2. * z1_pow4])  
        coef_2 = sum([4. * total_m_pow2 * zeta, 
                      4. * total_m * m_diff * pos_z1, 
                      -4. * m_diff * zeta * zeta_conj * pos_z1,
                      -2. * zeta * zeta_conj_pow2 * z1_pow2, 
                      4. * m_diff * z1_pow3, 2. * zeta * z1_pow4])
        coef_1 = sum([-8. * total_m * m_diff * zeta * pos_z1, 
                      -4. * m_diff_pow2 * z1_pow2, 
                      -4. * total_m_pow2 * z1_pow2, 
                      -4. * total_m * zeta * zeta_conj * z1_pow2, 
                      -4. * m_diff * zeta_conj * z1_pow3,
                      -zeta_conj_pow2 * z1_pow4, z1_pow3 * z1_pow3])
        coef_0 = sum([4. * m_diff_pow2 * zeta, 4. * total_m * m_diff * pos_z1, 
                      4. * m_diff * zeta * zeta_conj * pos_z1, 
                      2. * total_m * zeta_conj * z1_pow2, 
                      zeta * zeta_conj_pow2 * z1_pow2, 
                      -2. * m_diff * z1_pow3 - zeta * z1_pow4])
        coef_0 *= z1_pow2
        coefs_list = [coef_0, coef_1, coef_2, coef_3, coef_4, coef_5]
        return np.array(coefs_list).reshape(6)
        

    def _get_polynomial_roots_WM95(self, source_x, source_y):
        """roots of the polynomial"""
        polynomial_input = [self.mass_1, self.mass_2, self.separation, 
                            source_x, source_y]
        if polynomial_input != self._last_polynomial_input:
            polynomial = self._get_polynomial_WM95(source_x=source_x, 
                                                   source_y=source_y)
            self._polynomial_roots_WM95 = (
                    np.polynomial.polynomial.polyroots(polynomial))
            self._last_polynomial_input = polynomial_input
        return self._polynomial_roots_WM95

    def _polynomial_roots_ok_WM95(self, source_x, source_y):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots_WM95(source_x=source_x, 
                                                   source_y=source_y)
        component2 = self.mass_1 / np.conjugate(#Utils.complex_fsum(
#                                    list(roots) + [-self._position_z1_WM95]))
                     roots -self._position_z1_WM95)
        component3 = self.mass_2 / np.conjugate(#Utils.complex_fsum(
#                                    list(roots) + [-self._position_z2_WM95]))        
                roots -self._position_z2_WM95)
        #solutions = Utils.complex_fsum(
        #                            [self._zeta_WM95, component2, component3])
        solutions = self._zeta_WM95 + component2 + component3
        # This backs-up the lens equation.
        
        out = []
        for (i, root) in enumerate(roots):
            #min_distance_arg = np.argmin(abs(
            #                        Utils.complex_fsum([solutions, -root])**2))
            min_distance_arg = np.argmin(abs((solutions-root)**2))
            #print(min_distance_arg, abs((solutions-root)**2))
            if i == min_distance_arg:
                out.append(root)                 

        if len(out) not in [3, 5]:
            msg = 'CRITICAL ERROR - CONTACT CODE AUTHORS AND PROVIDE: {:} {:} {:} {:} {:}'
            txt = msg.format(repr(self.mass_1), repr(self.mass_2), 
                    repr(self.separation), repr(source_x), repr(source_y))
            # The repr() function gives absolute accuracy of float values allowing reproducing the results. 
            raise ValueError(txt)            
        return np.array(out)
        
    def _jacobian_determinant_ok_WM95(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._polynomial_roots_ok_WM95(
                                   source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        denominator_1 = self._position_z1_WM95 - roots_ok_bar
        add_1 = self.mass_1 / denominator_1**2
        denominator_2 = self._position_z2_WM95 - roots_ok_bar
        add_2 = self.mass_2 / denominator_2**2
        #derivative = Utils.complex_fsum([add_1, add_2])
        derivative = add_1 + add_2
        #print(derivative)
        #return Utils.complex_fsum([1., -derivative*np.conjugate(derivative)])
        return 1.-derivative*np.conjugate(derivative)

    def _signed_magnification_WM95(self, source_x, source_y):
        """signed magnification for each image separately"""
        return 1. / self._jacobian_determinant_ok_WM95(
                source_x=source_x, source_y=source_y)

    def _point_source_WM95(self, source_x, source_y):
        """calculate point source magnification using Witt & Mao 1995"""
        return fsum(abs(self._signed_magnification_WM95(
                source_x=source_x, source_y=source_y)))

    def _point_source_Witt_Mao_95(self, source_x, source_y):
        """calculate point source magnification using Witt & Mao 1995"""
        return self._point_source_WM95(
                source_x=source_x, source_y=source_y)

    def point_source_magnification(self, source_x, source_y):
        """calculates point source magnification for given position"""
        # Specify coordinate system convention!
        return self._point_source_Witt_Mao_95(
                source_x=source_x, source_y=source_y)

if __name__ == '__main__':
    s = 1.35
    q = 0.00578
    m1 = 1. / (1. + q)
    m2 = q / (1. + q)
    
    bl = BinaryLens(m1, m2, s)
    a = bl.point_source_magnification(1.38920106, 0.00189679)
    print(a)
    
