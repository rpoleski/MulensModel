import MulensModel
from astropy import units as u

l = MulensModel.Lens(n_components=2)
l.m_1 = 1.0 * u.M_sun
l.m_2 = 0.1 * u.M_jup
l.distance = 2. * u.kpc

s = MulensModel.Source()
s.distance = 8. * u.kpc

m = MulensModel.Model(lens=l, source=s)

# we have all what we need, hence we can access m.theta_E:
print(m.theta_E)

