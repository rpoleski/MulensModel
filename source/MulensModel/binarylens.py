#! /usr/bin/env python3

import numpy as np
from math import fsum, sqrt

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

        coef_5 = Utils.complex_fsum([z1_pow2, -zeta_conj_pow2])
        coef_4 = Utils.complex_fsum([-2. * total_m * zeta_conj, 
                       zeta * zeta_conj_pow2, -2. * m_diff * pos_z1, 
                       -zeta * z1_pow2])
        coef_3 = Utils.complex_fsum([4. * total_m * zeta * zeta_conj, 
                      4. * m_diff * zeta_conj * pos_z1, 
                      2. * zeta_conj_pow2 * z1_pow2, -2. * z1_pow4])  
        coef_2 = Utils.complex_fsum([4. * total_m_pow2 * zeta, 
                      4. * total_m * m_diff * pos_z1, 
                      -4. * m_diff * zeta * zeta_conj * pos_z1,
                      -2. * zeta * zeta_conj_pow2 * z1_pow2, 
                      4. * m_diff * z1_pow3, 2. * zeta * z1_pow4])
        coef_1 = Utils.complex_fsum([-8. * total_m * m_diff * zeta * pos_z1, 
                      -4. * m_diff_pow2 * z1_pow2, 
                      -4. * total_m_pow2 * z1_pow2, 
                      -4. * total_m * zeta * zeta_conj * z1_pow2, 
                      -4. * m_diff * zeta_conj * z1_pow3,
                      -zeta_conj_pow2 * z1_pow4, z1_pow3 * z1_pow3])
        coef_0 = Utils.complex_fsum([4. * m_diff_pow2 * zeta, 
                      4. * total_m * m_diff * pos_z1, 
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

    def _polynomial_roots_ok_WM95(self, source_x, source_y, return_distances=False):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots_WM95(source_x=source_x, 
                                                   source_y=source_y)
        component2 = self.mass_1 / np.conjugate(#Can we use Utils.complex_fsum()-like function here?
                                    roots -self._position_z1_WM95)
        component3 = self.mass_2 / np.conjugate(#Can we use Utils.complex_fsum()-like function here?
                                    roots -self._position_z2_WM95)
        solutions = self._zeta_WM95 + component2 + component3 #Can we use Utils.complex_fsum()-like function here?
        # This backs-up the lens equation.
        
        out = []
        distances = []
        for (i, root) in enumerate(roots):
            distances_from_root = abs((solutions-root)**2)
            min_distance_arg = np.argmin(distances_from_root) #Can we use Utils.complex_fsum()-like function here?
            if i == min_distance_arg:
                out.append(root)
                distances.append(distances_from_root[min_distance_arg])
        # The values in distances[] are a diagnostic on how good the numerical accuracy is.

        if len(out) not in [3, 5]:
            msg = 'CRITICAL ERROR - CONTACT CODE AUTHORS AND PROVIDE: {:} {:} {:} {:} {:}'
            txt = msg.format(repr(self.mass_1), repr(self.mass_2), 
                    repr(self.separation), repr(source_x), repr(source_y))
            # The repr() function gives absolute accuracy of float values allowing reproducing the results. 
            raise ValueError(txt)
        
        if return_distances:
            return (np.array(out), np.array(distances))
        else:
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
        derivative = add_1 + add_2 #Can we use Utils.complex_fsum()-like function here?
        return 1.-derivative*np.conjugate(derivative) #Can we use Utils.complex_fsum()-like function here?

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
        """calculates point source magnification for given position
        in a coordinate system where higher mass is at the origin
        and lower mass is at (separation, 0)"""
        return self._point_source_Witt_Mao_95(
                source_x=source_x + self.separation / 2., 
                source_y=source_y)

    def _get_magnification_w_plus(self, source_x, source_y, radius, 
                                  magnification_center=None):
        """Evaluates Gould (2008) eq. 7"""
        dx = [1., 0., -1., 0.]
        dy = [0., 1., 0., -1.]
        out = []
        for (i, dxval) in enumerate(dx):
            x = source_x + dxval * radius
            y = source_y + dy[i] * radius
            out.append(self.point_source_magnification(
                                              source_x=x, source_y=y))
        if magnification_center is None:
            magnification_center = self.point_source_magnification(
                                    source_x=source_x, source_y=source_y)
        return 0.25 * fsum(out) - magnification_center

    def _get_magnification_w_times(self, source_x, source_y, radius, 
                                  magnification_center=None):
        """Evaluates Gould (2008) eq. 8"""
        shift = radius / sqrt(2.)
        dx = [1., -1., -1., 1.]
        dy = [1., 1., -1., -1.]
        out = []
        for (i, dxval) in enumerate(dx):
            x = source_x + dxval * shift
            y = source_y + dy[i] * shift
            out.append(self.point_source_magnification(
                                              source_x=x, source_y=y))
        if magnification_center is None:
            magnification_center = self.point_source_magnification(
                                    source_x=source_x, source_y=source_y)
        return 0.25 * fsum(out) - magnification_center

    def hexadecapole_magnification(self, source_x, source_y, 
                                    rho, gamma, quadrupole=False):
        """hexadecpole approximation of binary-lens/finite-source 
        calculations - based on Gould 2008 ApJ 681, 1593"""
        # In this function, variables named a_* depict magnification.
        a_center = self.point_source_magnification(
                                    source_x=source_x, source_y=source_y)
        a_rho_half_plus = self._get_magnification_w_plus(source_x=source_x, 
                                            source_y=source_y, radius=0.5*rho, 
                                            magnification_center=a_center)
        a_rho_plus = self._get_magnification_w_plus(source_x=source_x, 
                                    source_y=source_y, radius=rho,
                                    magnification_center=a_center)
        
        # This is Gould 2008 eq. 9:
        a_2_rho_square = (16. * a_rho_half_plus - a_rho_plus) / 3. 
        
        # Gould 2008 eq. 6 (part 1/2):
        a_finite = a_center + a_2_rho_square * (1. - 0.2 * gamma)
        
        # At this point, a_finite is quadrupole approximation.
        if quadrupole:
            return a_finite
        
        a_rho_times = self._get_magnification_w_times(source_x=source_x, 
                                                source_y=source_y, radius=rho, 
                                                magnification_center=a_center)
        
        # This is Gould (2008) eq. 9:
        a_4_rho_power4 = 0.5 * (a_rho_plus + a_rho_times) - a_2_rho_square
        # This is Gould (2008) eq. 6 (part 2/2):
        a_finite += a_4_rho_power4 * (1. - 11. * gamma / 35.)
        return a_finite
        

if __name__ == '__main__':
    s = 1.35
    q = 0.00578
    m1 = 1. / (1. + q)
    m2 = q / (1. + q)
    x = 1.38920106 - s/2.
    y = 0.00189679
    rho = 0.001
    gamma = 0.5
    
    bl = BinaryLens(m1, m2, s)
    a = bl.point_source_magnification(x, y)
    print(a)
    aa = bl.hexadecapole_magnification(x, y, rho, gamma)
    print(aa)
    aaa = bl.hexadecapole_magnification(x, y, rho, gamma, quadrupole=True)
    print(aaa)
    