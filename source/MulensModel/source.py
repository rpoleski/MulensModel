from astropy import units as u

class LimbDarkeningCoeffs(object):
    def __init__(self, gamma_I=0.44, gamma_V=0.72, gamma_H=0.26):
        """
        Default values from Claret probably for a G dwarf.
        """
        self.gamma_I = gamma_I
        self.gamma_V = gamma_V
        self.gamma_H = gamma_H

default_ld_coeffs = LimbDarkeningCoeffs()

class Source(object):
    def __init__(distance=None, angular_size=None,
                 limb_darkening=None):
        self.distance=distance
        self.angular_size=angular_size
        if limb_darkening is None:
            self.limb_darkening = default_ld_coeffs

