from astropy import units as u

import MulensModel

#Point Lens
my_lens = MulensModel.Lens(mass=0.5*u.Mass, distance=6.e3*u.pc)
my_source = MulensModel.Source(distance=8.e3*u.pc)

point_lens = MulensModel.Model(lens=my_lens,source=my_source)
point_lens.plot_lightcurve() #plots magnification curve over +/-2*tE
point_lens.source.Imag = 18.
point_lens.blend.Imag = 20.
point_lens.plot_lightcurve(range=[-1.,1.])

#2 bodies - Version 1
l = MulensModel.Lens(n_components=2)
l.mass_1 = 1.0*u.solMass
l.mass_2 = 0.1*u.jupiterMass
l.distance = 2.e3*u.pc

s = MulensModel.Source()
s.distance = 8.e3*u.pc

m = MulensModel.Model(lens=l, source=s)

# we have all what we need, hence we can access m.theta_E:
print(m.theta_E)

#2 bodies - Version 2
two_body_lens = MulensModel.Lens(s=1.2,q=10.e-4)
two_body_lens.plot_caustics()
two_body_model = MulensModel.model(t0=5400.,u0=0.001,tE=14.*u.day,rho=0.001,
                                   s=two_body_lens.q, s=two_body_lens.s,
                                   alpha = 261.*u.deg)
two_body_model.plot_trajectory(caustics=True)

