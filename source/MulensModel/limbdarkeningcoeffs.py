from MulensModel.utils import Utils


class LimbDarkeningCoeffs(object):
    """
    Linear limb-darkening parameters. Both gamma and u conventions can be used. 
    The u convention is more frequently used in studies other in microlensing. 
    It has fixed flux at the center. An et al. 2002 (ApJ 572, 521,
    http://adsabs.harvard.edu/abs/2002ApJ...572..521A) introduced gamma convention:
    
    gamma = (2 * u) / (3 - u)
    
    u = 3 * gamma / (2 + gamma)
    
    Note that gamma convention has fixed total flux.
    """
    def __init__(self): #, gamma_I=0.44, gamma_V=0.72, gamma_H=0.26):
        self._gammas_for_band = dict()
        
    def set_limb_coef_gamma(self, bandpass, gamma):
        """set gamma coefficient for given band"""
        self._gammas_for_band[bandpass] = gamma
        
    def set_limb_coef_u(self, bandpass, u):
        """set u coefficient for given band"""
        self._gammas_for_band[bandpass] = Utils.u_to_gamma(u)
    
    def limb_coef_gamma(self, bandpass):
        """gives gamma coefficient for given band"""
        try:
            gamma = self._gammas_for_band[bandpass]
        except KeyError:
            msg = ('No limb darkening coefficient for bandpass {:}. Most ' + 
                   'probably you have set the filter for a dataset but have' +
                   ' not set LD coef for this filter')
            raise KeyError(msg.format(bandpass))
        return gamma
    
    def limb_coef_u(self, bandpass):
        """gives u coefficient for given band"""
        gamma = self.limb_coef_gamma(bandpass=bandpass)
        return Utils.gamma_to_u(gamma)
        
    def weighted_limb_coef_gamma(self, weights):
        """get weighted limb darkening coefficient
        weights is a distionary e.g. weights = {'I': 1.5, 'V': 1.} """
        if not isinstance(weights, dict):
            raise TypeError("LimbDarkeningCoeffs.weighted_limb_coef_gamma() " +
                    "parameter has to be dict, not {:}".format(type(weights)))
        gamma_sum = 0.
        weight_sum = 0.
        for (band, weight) in weights.items():
            try:
                gamma_sum += self.limb_coef_gamma(band)
            except KeyError:
                msg = "The bandpass {:} was not set for limb darkening"
                raise KeyError(msg.format(band))
            weight_sum += weight
        return gamma_sum / weight_sum
        
