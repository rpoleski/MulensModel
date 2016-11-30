from astropy import units as u

import MulensModel

#Point Lens: Implementation of MulensModel
my_lens = MulensModel.Lens(mass=0.5*u.Mass, distance=6.e3*u.pc)
my_source = MulensModel.Source(distance=8.e3*u.pc)

point_lens = MulensModel.Model(u_0=0.1, lens=my_lens, source=my_source)
point_lens.plot_lightcurve() #plots magnification curve over +/-2*tE
point_lens.source.I_mag = 18.
point_lens.blend.I_mag = 20. # what type of object is MulensModel.Model.blend ???
point_lens.plot_lightcurve(range=[-1.,1.])

#2 bodies - Version 1: Implementation of MulensSystem
two_body_lens = MulensModel.Lens()
two_body_lens.mass_1 = 1.0 * u.solMass
two_body_lens.mass_2 = 0.1 * u.jupiterMass
two_body_lens.distance = 2.e3 * u.pc

source = MulensModel.Source()
source.distance = 8.e3*u.pc

two_body_system = MulensModel.MulensSystem(lens=two_body_lens, source=source)
print(two_body_system.theta_E)

#2 bodies - Version 2: Implementation of MulensModel
two_body_lens_v2 = MulensModel.Lens(s=1.2, q=10.e-4)
two_body_lens_v2.plot_caustics()
two_body_model = MulensModel.Model(
    t_0=5400., u_0=0.001, t_E=14.*u.day, rho=0.001, s=two_body_lens_v2.q, 
    s=two_body_lens_v2.s, alpha=261.*u.deg)
two_body_model.plot_trajectory(caustics=True)

