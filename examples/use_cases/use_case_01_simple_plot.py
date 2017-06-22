"""
Use Case 01: Define and plot a 2-body microlensing model
"""
import matplotlib.pyplot as pl
from astropy import units as u

import MulensModel

model = MulensModel.Model()

print(model.parameters)

model.set_parameters(
    t_0=2457603.1, u_0=0.23, t_E=45*u.day, rho=0.001, alpha=130.23*u.deg, s=1.3, 
    q=0.3)
print(model.parameters)

pl.figure()
pl.subplot(2, 1, 1)
pl.title('Magnification Curve')
model.plot_magnification()

pl.subplot(2, 1, 2)
pl.title('Caustic Structure & Trajectory')
model.plot_trajectory(caustics=True)

pl.show()

