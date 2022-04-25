from MulensModel.utils import Utils


class LimbDarkeningCoeffs(object):
    """
    Linear limb-darkening parameters. Both *gamma* and *u* conventions
    can be used.  The *u* convention is more frequently used in
    studies other than microlensing.  It has fixed flux at the center.
    `An et al. 2002 (ApJ 572, 521)
    <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_
    introduced the *gamma* convention:

    gamma = (2 * u) / (3 - u)

    u = 3 * gamma / (2 + gamma)

    Note that the gamma convention has fixed total flux.

    """

    def __init__(self):
        self._gammas_for_band = dict()

    def set_limb_coeff_gamma(self, bandpass, gamma):
        """
        Remembers limb darkening gamma coefficient for given band.

        Parameters :
            bandpass: *str*
                Name of the filter.

            gamma: *float*
                The value of gamma coefficient.

        """
        self._gammas_for_band[bandpass] = gamma

    def get_limb_coeff_gamma(self, bandpass):
        """
        Gives limb darkening gamma coefficient for given band.

        Parameters :
            bandpass: *str*
                Name of the filter.

        Returns :
            gamma: *float*
                The value of gamma coefficient.

        """
        try:
            gamma = self._gammas_for_band[bandpass]
        except KeyError:
            msg = ('No limb darkening coefficient for bandpass {:}. Most ' +
                   'probably you have set the filter for a dataset but have' +
                   ' not set LD coeff for this filter')
            raise KeyError(msg.format(bandpass))
        return gamma

    def set_limb_coeff_u(self, bandpass, u):
        """
        Remembers limb darkening *u* coefficient for given band

        Parameters :
            bandpass: *str*
                Name of the filter.

            u: *float*
                The value of *u* coefficient.

        """
        self._gammas_for_band[bandpass] = Utils.u_to_gamma(u)

    def get_limb_coeff_u(self, bandpass):
        """
        Gives limb darkening *u* coefficient for given band.

        Parameters :
            bandpass: *str*
                Name of the filter.

        Returns :
            u: *float*
                The value of *u* coefficient.

        """
        gamma = self.get_limb_coeff_gamma(bandpass=bandpass)
        return Utils.gamma_to_u(gamma)

    def get_weighted_limb_coeff_gamma(self, weights):
        """
        Get weighted limb darkening coefficient in gamma space.

        Parameters :
            weights: *dict*
                A dictionary that for every band (keys; *str* type)
                gives its relative weight (value; *float* type), e.g.,
                ``weights = {'I': 1.5, 'V': 1.}`` will return gamma
                coefficient in the case when *I* band contributes 1.5
                more than *V* band. Note that for each band used you
                have to first set to coefficient.

        Returns :
            gamma: *float*
                The value of weighted gamma.

        """

        if not isinstance(weights, dict):
            raise TypeError(
                "LimbDarkeningCoeffs.get_weighted_limb_coeff_gamma() " +
                "parameter has to be dict, not {:}".format(type(weights)))
        gamma_sum = 0.
        weight_sum = 0.
        for (band, weight) in weights.items():
            try:
                gamma_sum += self.get_limb_coeff_gamma(band)
            except KeyError:
                msg = "The bandpass {:} was not set for limb darkening"
                raise KeyError(msg.format(band))
            weight_sum += weight
        return gamma_sum / weight_sum
