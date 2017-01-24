from astropy import units as u


class LimbDarkeningCoeffs(object):
    """
    Linear limb-darkening parameters. Defined as gamma = (2 * u) / (3 - u)
    """
    # Possibly this needs to be a dictionary or just more
    # extensible. Think MOA limb-darkening.
    def __init__(self, gamma_I=0.44, gamma_V=0.72, gamma_H=0.26):
        """
        Set the limb-darkening coefficients.

        Default values from Claret probably for a G dwarf.
        """
        self.gamma_I = gamma_I
        self.gamma_V = gamma_V
        self.gamma_H = gamma_H

