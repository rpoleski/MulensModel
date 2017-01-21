import numpy as np
from math import fsum


class BinaryLensEquation(object):
    """binary lens equation - its solutions, images, parities etc.
    Calculations are based on Witt & Mao (1995).
    """
    def __init__(self, mass_1=None, mass_2=None, separation=None, source_x=None, source_y=None):
        """mass_1, mass_2, and separation are relative to some mass (and corresponding Einstein radius)"""
        self.mass_1 = mass_1
        self.mass_2 = mass_2
        self.separation = separation
        self.source_x = source_x
        self.source_y = source_y
        self._total_mass = None
        self._mass_difference = None
        self._position_z1 = None
        self._position_z2 = None
        self._zeta = None
        self._polynomial_coefs = None
        self._polynomial_roots = None
        self._polynomial_roots_verified = None
        self._jacobian_determinant_verified = None
        self._signed_magnification = None
        
    def _calculate_variables(self):
        """calculates values of variables needed for polynomial coefficients"""
        self._total_mass = 0.5 * (self.mass_1 + self.mass_2)
        self._mass_difference = 0.5 * (self.mass_2 - self.mass_1)
        self._position_z1 = -0.5 * self.separation + 0.j
        self._position_z2 = 0.5 * self.separation + 0.j
        self._zeta = self.source_x - self.separation + self.source_y * 1.j
        
    def _get_polynomial_coefs(self):
        """calculate coefficients of the polynomial"""
        if None in [self._total_mass, self._mass_difference, self._position_z1, self._position_z2, self._zeta]:
            self._calculate_variables()
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m
        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff
        pos_z1 = self._position_z1
        pos_z2 = self._position_z2
        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2
        zeta = self._zeta
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
        self._polynomial_coefs = np.array(coefs_list).reshape(6)

    @property
    def polynomial_coefs(self):
        """coefficients for lens polynomial - convention as in numpy.polynomial.polynomial"""
        if self._polynomial_coefs is None:
            self._get_polynomial_coefs()
        return self._polynomial_coefs

    @property 
    def polynomial_roots(self):
        """roots of the polynomial"""
        if self._polynomial_roots is None:
            self._polynomial_roots = np.polynomial.polynomial.polyroots(self.polynomial_coefs)
        return self._polynomial_roots
    
    @property
    def polynomial_roots_verified(self):
        """roots of polynomial that are also roots of the lens equation"""
        if self._polynomial_roots_verified is None:
            component2 = self.mass_1 / np.conjugate(sum([self.polynomial_roots, -self._position_z1]))
            component3 = self.mass_2 / np.conjugate(sum([self.polynomial_roots, -self._position_z2]))
            solutions = sum([self._zeta, component2, component3]) # This backs-up the lens equation
            out = []
            for i in range(len(self.polynomial_roots)):
                min_distance_arg = np.argmin(abs((solutions-self.polynomial_roots[i])**2))
#                print(i, min_distance_arg, abs((solutions-self.polynomial_roots[i])**2)[min_distance_arg])
                if i == min_distance_arg:
                    out.append(self.polynomial_roots[i])
            self._polynomial_roots_verified = np.array(out)
            if len(self._polynomial_roots_verified) not in [3, 5]:
                msg = 'CRITICAL ERROR - CONTACT CODE AUTHORS AND PROVIDE: {:} {:} {:} {:} {:}'
                txt = msg.format(repr(self.mass_1), repr(self.mass_2), repr(self.separation), repr(self.source_x.value[0]), repr(self.source_y.value[0]))
                # The repr() function gives absolute accuracy of float values allowing reproducing the results. 
                raise ValueError(txt)
        return self._polynomial_roots_verified
        
    @property
    def jacobian_determinant_verified(self):
        """determinants of lens equation Jacobian"""
        if self._jacobian_determinant_verified is None:
            roots_ok_conj = np.conjugate(self.polynomial_roots_verified)
            denominator_1 = self._position_z1 - roots_ok_conj
            denominator_1 *= denominator_1
            denominator_2 = self._position_z2 - roots_ok_conj
            denominator_2 *= denominator_2
            derivative = sum([self.mass_1 / denominator_1, self.mass_2 / denominator_2])
            self._jacobian_determinant_verified = sum([1., -derivative * np.conjugate(derivative)])
        return self._jacobian_determinant_verified
        
    @property
    def signed_magnification(self):
        """signed magnification for each image separately"""
        if self._signed_magnification is None:
            self._signed_magnification = 1. / self.jacobian_determinant_verified
        return self._signed_magnification
         
    @property
    def total_magnification(self):
        """total magnification of all images"""
        return [sum(abs(self.signed_magnification))]

