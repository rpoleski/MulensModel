import numpy as np

import MulensModel


t = np.arange(7000., 7600., 0.1)
q = np.zeros(t.size)
e = np.zeros(t.size)

d = MulensModel.MulensData(np.array(t, q, e), satellite='K2')

print(d.satellite)

m = MulensModel.Model()
m.parallax(satellite=True, earth_orbital=False, topocentric=False) 
#default is all True



