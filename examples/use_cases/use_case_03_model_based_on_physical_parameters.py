"""
Use Case 03: Define a Model Based on a Physical System

"""
from astropy import units as u
import matplotlib.pyplot as pl

import MulensModel

#Define a Point Lens Model via MulensSystem
my_lens = MulensModel.Lens(mass=0.5*u.Mass, distance=6.e3*u.pc)
my_source = MulensModel.Source(distance=8.e3*u.pc)
my_system = MulensModel.MulensSystem(lens=my_lens, source=my_source)
print(my_system)

pl.figure()
pl.title('Magnification Curve')
my_system.plot_magnification(u0=0.3)

my_system.mu_rel = 3. * u.mas / u.yr
my_model = MulensModel.Model(t0=7620., u0=0.3, tE=my_system.tE)

#Define a 2-body model Version 1: Implementation via MulensSystem
two_body_lens = MulensModel.Lens()
two_body_lens.mass_1 = 1.0 * u.solMass
two_body_lens.mass_2 = 0.1 * u.jupiterMass
two_body_lens.distance = 2.e3 * u.pc

source = MulensModel.Source()
source.distance = 8.e3*u.pc

two_body_system = MulensModel.MulensSystem(lens=two_body_lens, source=source)
print(two_body_system)

#Define a 2-body model Version 2: Implementation via MulensModel
two_body_lens_v2 = MulensModel.Lens(s=1.2, q=10.e-4)
pl.figure()
pl.title('2-body caustic structure')
two_body_lens_v2.plot_caustics()

two_body_model = MulensModel.Model(
    t_0=5400., u_0=0.001, t_E=14.*u.day, rho=0.001, s=two_body_lens_v2.q, 
    s=two_body_lens_v2.s, alpha=261.*u.deg)
pl.figure()
pl.title('Caustics & Trajectory')
two_body_model.plot_trajectory(caustics=True)

pl.show()
