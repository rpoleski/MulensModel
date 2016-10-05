from astropy import units as u

import MulensModel


my_lens = MulensModel.Lens(mass=0.5*u.Mass, distance=6.e3*u.pc)
my_source = MulensModel.Source(distance=8.e3*u.pc)

point_lens = MulensModel.Model(lens=my_lens,source=my_source)
point_lens.plot_lightcurve() #plots magnification curve over +/-2*tE
point_lens.source.Imag = 18.
point_lens.blend.Imag = 20.
point_lens.plot_lightcurve(range=[-1.,1.])


l = MulensModel.Lens(n_components=2)
l.mass_1 = 1.0*u.solMass
l.mass_2 = 0.1*u.jupiterMass
l.distance = 2.e3*u.pc

s = MulensModel.Source()
s.distance = 8.e3*u.pc

m = MulensModel.Model(lens=l, source=s)

# we have all what we need, hence we can access m.theta_E:
print(m.theta_E)

