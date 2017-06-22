from astropy import units as u

from MulensModel.model import Model
from MulensModel.mulensparallaxvector import MulensParallaxVector

"""
Use Case 13: Set the parallax parameters of a model.
"""

model = Model()

#Event Coordinates
model.ra = 270. * u.deg
model.dec = -28.5 * u.deg

#Parallax parameters
model.parallax(earth_orbital=True)
model.pi_E = MulensParallaxVector(0.2, 0.4, ref="NorthEast") 
    #NorthEast is default, the other choice is "par_perp"
model.t_0_par = 2457000.

print(model.pi_E)



