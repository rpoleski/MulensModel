import matplotlib.pyplot as pl
from astropy import units as u

import MulensModel


model = MulensModel.Model()

model.parameters(
    t_0=7603.1, u_0=0.23, t_E=45*u.day, rho=0.001, alpha=130.23*u.deg, s=1.3, 
    q=0.3)

pl.subplot(2, 1, 1)
pl.plot(model.time, model.magnification)

pl.subplot(2, 1, 2)
pl.plot(model.caustic.x, model.caustic.y)
pl.plot(model.trajectory.x, model.trajectory.y)

pl.show()

