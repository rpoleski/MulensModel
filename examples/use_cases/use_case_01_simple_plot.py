import MulensModel
import matplotlib.pyplot as pl
from astropy import units as u

m=MulensModel.Model()

m.parameters(t_0=7603.1, u_0=0.23, t_E=45*u.day, rho=0.001, alpha=130.23*u.deg, s=1.3, q=0.3)

pl.subplot(2, 1, 1)
pl.plot(m.t, m.mag)

pl.subplot(2, 1, 2)
pl.plot(m.caustic.x, m.caustic.y)
pl.plot(m.trajectory.x, m.trajectory.y)

pl.show()

