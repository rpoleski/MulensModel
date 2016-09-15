from astropy import units as u

import MulensModel

m = MulensModel.Model()

m.ra = 270. * u.deg
m.dec = -28.5 * u.deg
m.parallax(True)
m.pi_E((0.2, 0.4), ref="NorthEast") # this is default, the other choice is "par_perp"

m.t_0_par = 7000. # is it OK, or we want here instance of Time class from astropy ???

print(m.pi_E)



