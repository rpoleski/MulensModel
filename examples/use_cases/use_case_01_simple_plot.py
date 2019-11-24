"""
Use Case 01: Define and plot a 2-body microlensing model
"""
import matplotlib.pyplot as plt
from astropy import units as u

import MulensModel as mm


model = mm.Model(
    {'t_0': 2457603.1, 'u_0': 0.23, 't_E': 45*u.day,
     'alpha': 130.23*u.deg, 's': 1.3, 'q': 0.3})
print(model.parameters)

plt.figure()
plt.subplot(2, 1, 1)
plt.title('Magnification Curve')
model.plot_magnification()

plt.subplot(2, 1, 2)
plt.title('Caustic Structure & Trajectory')
model.plot_trajectory(caustics=True)

plt.show()
