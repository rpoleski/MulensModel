"""
Use Case 13: Set the parallax parameters of a model.
"""
from astropy import units as u

import MulensModel as mm


model = mm.Model(
    {'t_0': 2457005., 'u_0': 0.1, 't_E': 30., 'pi_E_N': 0.2, 'pi_E_E': 0.4, 't_0_par': 2457000.},
    ra=270.*u.deg, dec=-28.5*u.deg)
model.parallax(earth_orbital=True)

print(model.parameters.pi_E_N)
