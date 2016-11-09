from astropy import units as u

import MulensModel

model = MulensModel.Model()

model.ra = 270. * u.deg
model.dec = -28.5 * u.deg
model.parallax(True)
model.pi_E((0.2, 0.4), ref="NorthEast") # this is default, the other choice is "par_perp"

model.t_0_par = 7000. # is it OK, or we want here instance of Time class from astropy ???

print(model.pi_E)



