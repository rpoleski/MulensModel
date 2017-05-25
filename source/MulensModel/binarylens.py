#! /usr/bin/env python3

import sys
import numpy as np
from math import fsum, sqrt

from MulensModel.utils import Utils
# for VBBL imports see self._vbbl_imported below


class BinaryLens(object):
    """
    binary lens equation - its solutions, images, parities, magnifications, etc.
    
    The binary lens equation is a 5th order complex polynomial.
    """
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
        self._vbbl_imported = False

    def _calculate_variables(self, source_x, source_y):
        """calculates values of constants needed for polynomial coefficients"""
        self._total_mass = 0.5 * (self.mass_1 + self.mass_2) # This is total_mass in WM95 paper.
        self._mass_difference = 0.5 * (self.mass_2 - self.mass_1)
        self._position_z1_WM95 = -0.5 * self.separation + 0.j
        self._position_z2_WM95 = 0.5 * self.separation + 0.j
        self._zeta_WM95 = source_x - self.separation + source_y * 1.j

    def _get_polynomial_WM95(self, source_x, source_y):
        """calculate coefficients of the polynomial"""
        #Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff

        pos_z1 = self._position_z1_WM95

        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2

        zeta = self._zeta_WM95
        zeta_conj = np.conjugate(zeta)
        zeta_conj_pow2 = zeta_conj * zeta_conj

        #Calculate the coefficients of the 5th order complex polynomial
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

        #Return the coefficients of the polynomial
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

        #If the lens equation is solved correctly, there should be
        #either 3 or 5 solutions (corresponding to 3 or 5 images)
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
        """
        calculates point source magnification for given position
        in a coordinate system where the center of mass is at origin 
        and both masses are on X axis with higher mass at negative X;
        this means that the higher mass is at (X, Y)=(-s*q/(1+q), 0) and
        the lower mass is at (s/(1+q), 0)
        """
        x_shift = self.separation * (0.5 + 
                                    self.mass_2 / (self.mass_1 + self.mass_2))
        # We need to add this because WM95 use geometric center as an origin 
        # of their coordinate system.
        return self._point_source_Witt_Mao_95(
                #source_x=source_x + self.separation / 2. + self.separation * self.mass_2 / (self.mass_1 + self.mass_2), 
                source_x=source_x+x_shift, 
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

    def hexadecapole_magnification(self, source_x, source_y, rho, gamma, 
                                  quadrupole=False, all_approximations=False):
        """hexadecpole approximation of binary-lens/finite-source 
        calculations - based on Gould 2008 ApJ 681, 1593
        
        for coordinate system convention see point_source_magnification()"""
        # In this function, variables named a_* depict magnification.
        if quadrupole and all_approximations:
            msg = 'Inconsisient parameters of {:}'
            raise ValueError(msg.format('BinaryLens.hexadecapole_magnification()'))
        
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
        a_quadrupole = a_center + a_2_rho_square * (1. - 0.2 * gamma)
        
        # At this point is quadrupole approximation is finished
        if quadrupole:
            return a_quadrupole
        
        a_rho_times = self._get_magnification_w_times(source_x=source_x, 
                                                source_y=source_y, radius=rho, 
                                                magnification_center=a_center)
        
        # This is Gould (2008) eq. 9:
        a_4_rho_power4 = 0.5 * (a_rho_plus + a_rho_times) - a_2_rho_square
        # This is Gould (2008) eq. 6 (part 2/2):
        a_add = a_4_rho_power4 * (1. - 11. * gamma / 35.)
        a_hexadecapole = a_quadrupole + a_add
        
        if all_approximations:
            return (a_hexadecapole, a_quadrupole, a_center)
        else:
            return a_hexadecapole
            
    def vbbl_magnification(self, source_x, source_y, rho, gamma=0., accuracy=0.001):
        """calculate magnification using VBBL code that implements advanced 
        contour integration algorithm presented by
        Bozza 2010 MNRAS, 408, 2188
        http://adsabs.harvard.edu/abs/2010MNRAS.408.2188B
        http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm
        """
        if not self._vbbl_imported:
            try:
                from VBBL import VBBinaryLensingLibrary
            except ImportError:
                print("VBBL library could not be imported")
                print("You must compile VBBL library and add it to PYTHONPATH")
                print("Exit...")
                sys.exit(1)
            self._vbbl_imported = True
            self._vbbllib = VBBinaryLensingLibrary.VBBinaryLensing()
        
        s = float(self.separation)
        q = float(self.mass_2 / self.mass_1)
        x = float(source_x)
        y = float(source_y)
        rho = float(rho)
        if gamma != 0.:
            raise ValueError("The case of gamma != 0 is not yet coded.")
            # Note that VBBL uses u coefficient, not Gamma.
        gamma = float(gamma)
        accuracy = float(accuracy) # Note that this accuracy is not guaranteed.
        assert accuracy > 0., "VBBL requires accuracy > 0 e.g. 0.01 or 0.001;\n{:} was provided".format(accuracy)
        
        return self._vbbllib.BinaryMagDark(s, q, x, y, rho, gamma, accuracy)
        # To get the image positions from VBBL, following C++ code has to be run:
        #  _sols *Images;
        #  Mag=VBBL.BinaryMag(s, q, y1, y2, rho, accuracy, &Images);
        #  for(_curve *c=Images->first; c; c=c->next) {
        #      for(_point *p=c->first; p; p=p->next) {
        #          // Here p->x1 and p->x2 give next point 
        #      }
        #  }
        #  delete Images;
