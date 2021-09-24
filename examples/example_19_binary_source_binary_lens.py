"""
Example usage of the binary-lens binary-source model.
"""
import matplotlib.pyplot as plt
import os

import MulensModel as mm


t_0_1 = 5370.
t_0_2 = 5390.
u_0_1 = -0.15
u_0_2 = 0.2
rho_1 = 0.01
t_star_2 = 0.45

t_start = 5365.
t_stop = 5398.

common = {'t_E': 25., 's': 1.1, 'q': 0.1, 'alpha': 10.}
parameters = {'t_0_1': t_0_1, 't_0_2': t_0_2, 'u_0_1': u_0_1, 'u_0_2': u_0_2,
              'rho_1': rho_1, 't_star_2': t_star_2, **common}
parameters_1 = {'t_0': t_0_1, 'u_0': u_0_1, 'rho': rho_1, **common}
parameters_2 = {'t_0': t_0_2, 'u_0': u_0_2, 't_star': t_star_2, **common}

model = mm.Model(parameters)
model_1 = mm.Model(parameters_1)
model_2 = mm.Model(parameters_2)

model.set_magnification_methods([t_start, 'VBBL', t_stop])
model_1.set_magnification_methods([t_start, 'VBBL', t_stop])
model_2.set_magnification_methods([t_start, 'VBBL', t_stop])

# Make magnification plots:
plot_kwargs = {'t_start': t_start, 't_stop': t_stop}

model.plot_magnification(source_flux_ratio=1., label='ratio = 1',
                         **plot_kwargs)
model.plot_magnification(source_flux_ratio=5., label='ratio = 5',
                         **plot_kwargs)
model.plot_magnification(source_flux_ratio=0.2, label='ratio = 0.2',
                         **plot_kwargs)
model_1.plot_magnification(label='source 1', ls='--', **plot_kwargs)
model_2.plot_magnification(label='source 2', ls='--', **plot_kwargs)

plt.legend()
plt.show()

# Make trajectory plots:
model.plot_trajectory(caustics=True, label='both_sources', **plot_kwargs)
model_2.plot_trajectory(label='source 2', ls='--', lw=4, **plot_kwargs)
model_1.plot_trajectory(label='source 1', ls='--', lw=4, **plot_kwargs)

plt.legend()
plt.show()

#
file_name = os.path.join(mm.DATA_PATH, "photometry_files", "OB08092",
                         "phot_ob08092_O4.dat")
data = mm.MulensData(file_name=file_name)
