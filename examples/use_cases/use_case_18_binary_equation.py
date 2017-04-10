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

time_1 = 2455746.
time_2 = 2455746.6
time_3 = 2455746.7
time_4 = 2455747.
time_5 = 2455747.15
time_6 = 2455748.

model.set_magnification_method(
    {'Quadrupole':[(time_1, time_2), (time_5, time_6)],
     'Hexadecapole':[(time_2, time_3), (time_4, time_5)], 
     'VBBL':(time_3, time_4)})
# We would also like to specify a single method that is used for every
# epoch without specifying begin and end. This would be useful when we
# have e.g. many datasets and don't want to check min and max for
# every one. --> Write a use case.


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
