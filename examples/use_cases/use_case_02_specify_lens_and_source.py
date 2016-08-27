import MulensModel
from astropy import units as u

l = MulensModel.Lens(n_components=2)
l.s = 1.3
l.q = 0.1

s = MulensModel.Source()
s.trajectory(t_0=7263.16, u_0=0.2, alpha=130.12*u.deg)

# now combine lens, source and add remaining parameters
m = MulensModel.Model(lens=l, source=s)
m.t_E = 25. * u.day
m.rho = 0.001

