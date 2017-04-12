"""
Use Case 18: Binary Magnification Equation

Change how the binary magnification is calculated based on 
- the time stamp of the model
- the magnificaiton of the model
"""
from astropy import units as u
import numpy as np
import matplotlib.pyplot as pl

from MulensModel import Model

#Initialize the model 
t_0 = 2455747.049357
t_E = 21.6796
model = Model(
    t_0=t_0, u_0=0.00352, t_E=t_E, rho=0.001647, alpha=41.35*u.deg, s=0.5486, 
    q=0.00532)

#times to calculate the magnification
times = np.arange(t_0 - 1., t_0 + 1., 0.001)

#Calculate the magnification using different magnification calculations
default_magnification = model.magnification(times)

times = np.array([46., 46.6, 46.7, 47., 47.15, 48.]) + 2455700. # There are 6 of them.
methods = ['Quadrupole', 'Hexadecapole', 'VBBL', 'Hexadecapole', 'Quadrupole'] # There are 5 of them.

model.set_magnification_methods(dividing_epochs=times, methods=methods)

#NOT IMPLEMENTED: Set times of caustic crossings. Use different
#magnification calculation based on number of source radii from the
#crossing.

accurate_magnification_2 = model.magnification(times)

#Plot the differences
pl.figure()
pl.title('Magnification Curves')
pl.plot(times, default_magnification, color='black', label='Default Magnification')
pl.plot(times, accurate_magnification, color='green', label='User Specified Equations')
pl.xlabel('Times')
pl.ylabel('Magnification')
pl.legend(loc='best')

pl.figure()
pl.title('Difference in Magnification Curves')
pl.plot(times, default_magnification - accurate_magnification, color='green')
pl.xlabel('Times')
pl.ylabel('Magnification Difference')
pl.legend(loc='best')

pl.show()
