import numpy as np

class ModelIMF(object):
    """
    Class for calculating IMF
    """
    def __init__(self):
        self._alpha = np.array([0.3, 1.3, 2.3])
        #self._mass_boundary = np.array([0.01, 0.08, 0.5, 150.0])
        # Let's shift it to allow less massive lenses for K2C9-1:
        self._mass_boundary = np.array([0.001, 0.08, 0.5, 150.0])

        delta_alpha = self._alpha[1:] - self._alpha[:-1]
        boundary = self._mass_boundary[1:-1]
        self._factors = np.insert(np.cumprod(boundary**delta_alpha), 0, 1.)

    def get_relative_probability(self, value):
        """XXX"""
        if value < self._mass_boundary[0] or value > self._mass_boundary[-1]:
            return -np.inf
        x = np.digitize(value, self._mass_boundary) - 1
        return self._factors[x] * value**-self._alpha[x]

    def get_ln_relative_probability(self, value):
        """XXX"""
        if value < self._mass_boundary[0] or value > self._mass_boundary[-1]:
            return -np.inf
        x = np.digitize(value, self._mass_boundary) - 1
        return np.log(self._factors[x]) - self._alpha[x] * np.log(value)

if __name__ == '__main__':
    m = ModelIMF()
    for i in np.arange(0.015, 0.9, 0.01):
        print(i, m.get_relative_probability(i))
