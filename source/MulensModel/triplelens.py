import numpy as np


class TripleLens(object):
    """
    The triple lens equation - its solutions, images, parities,
    magnifications, etc.

    Arguments :
        mass_1: *float*
            Mass of the primary as a fraction of the total mass.

        mass_2: *float*
            Mass of the secondary as a fraction of the total mass.

        mass_3: *float*
            Mass of the tertiary as a fraction of the total mass.
            XXX NOTE - add conventions

        separation_21: *float*
            Separation between bodies 2 and 1 as a fraction of
            the Einstein ring of the whole system.

        separation_31: *float*
            Separation between bodies 2 and 1 as a fraction of
            the Einstein ring of the whole system.

        psi: *float*
            XXX NOTE - description

    Note: masses 1, 2, and 3 may be defined as a fraction of some other
    mass than the total mass. This is possible but not recommended -
    make sure you know what you're doing before you start using this
    possibility.
    """
    def __init__(self, mass_1=None, mass_2=None, mass_3=None,
                 separation_21=None, separation_31=None, psi=None):
        self._mass_1 = mass_1
        self._mass_2 = mass_2
        self._mass_3 = mass_3
        self._separation_21 = separation_21
        self._separation_31 = separation_31
        self._psi = psi

    def get_point_source_magnification(self, source_x, source_y):
        """
        XXX

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

        Returns :
            magnification: *float*
                Point source magnification.
        """
        pass
        # use https://arxiv.org/abs/astro-ph/0202294

# XXX :
#    def get_hexadecapole_magnification(self, source_x, source_y, rho, gamma,
#                               quadrupole=False, all_approximations=False):
