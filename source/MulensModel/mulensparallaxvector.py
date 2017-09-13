import numpy as np

class MulensParallaxVector(object):
    """
    A class for the microlens parallax, which is a vector. May be
    specified either relative to the sky ("NorthEast") or relative to
    the binary lens axis ("ParPerp"). "NorthEast" is default.
    """
    def __init__(self, pi_E_1=None, pi_E_2=None, pi_E=None, ref=None):
        """
        Specify the components of the parallax vector. Either by
        specifying the components individually (i.e. as pi_E_1 and
        pi_E_2) or by specifying pi_E as a list, tuple, np.ndarray
        with 2 elements.

        May be specified either relative to the sky (ref="NorthEast")
        or relative to the binary lens axis (ref="ParPerp"). "NorthEast" is
        default.
        """
        #Set reference frame
        if ref is None:
            self.ref = "NorthEast"
        else:
            self.ref = ref

        #Store the parallax vector
        if pi_E_1 is not None:
            if pi_E_2 is not None:
                self.vector = np.array((pi_E_1, pi_E_2))
            else:
                raise AttributeError('pi_E has 2 components')
        else:
            if pi_E is not None:
                pi_E = np.array(pi_E)
                if pi_E.size is 2:
                    self.vector = pi_E
                else:
                    raise AttributeError('pi_E has 2 components')

    def __repr__(self):
        """A nice way to print the MulensParallaxVector"""

        parallax_str = "pi_E"
        if self.ref == "NorthEast":
            parallax_str= "(pi_E_N, pi_E_E)"
        if self.ref == "ParPerp":
            parallax_str = "(pi_E_par, pi_E_perp)"

        return "{0} = ({1}, {2})".format(parallax_str, self.vector[0], 
                                         self.vector[1])
