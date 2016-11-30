from astropy import units as u

import MulensModel

point_lens_2 = MulensModel.Model(t_0=7600., u_0=0.1, t_E=25.*u.day)
print(point_lens_2.parameters())
point_lens_2.plot_lightcurve()

planet_model = MulensModel.Model(
    t_0=7650., u_0=0.3, t_E=32., rho=0.0005, s=0.8, q=10.e-3, alpha=32.*u.deg)
print(planet_model.lens)
planet_model.source = MulensModel.Source(distance=8.e3*u.pc,
                                         angular_size=0.001*u.mas)
print(planet_model.thetaE)
#Note there is probably a way to get microarcseconds in astropy.units.cds
