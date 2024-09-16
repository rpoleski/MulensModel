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

common = {'t_E': 25., 's': 1.1, 'q': 0.1, 'alpha': 190.}
parameters = {'t_0_1': t_0_1, 't_0_2': t_0_2, 'u_0_1': u_0_1, 'u_0_2': u_0_2,
              'rho_1': rho_1, 't_star_2': t_star_2, **common}
parameters_1 = {'t_0': t_0_1, 'u_0': u_0_1, 'rho': rho_1, **common}
parameters_2 = {'t_0': t_0_2, 'u_0': u_0_2, 't_star': t_star_2, **common}

# plotting ranges:
t_start = 5362.
t_stop = 5400.

model = mm.Model(parameters)
model_1 = mm.Model(parameters_1)
model_2 = mm.Model(parameters_2)

model.set_magnification_methods([t_start, 'VBBL', t_stop])
model_1.set_magnification_methods([t_start, 'VBBL', t_stop])
model_2.set_magnification_methods([t_start, 'VBBL', t_stop])

# Make magnification plots:
plt.figure()
plot_kwargs = {'t_start': t_start, 't_stop': t_stop}

model.plot_magnification(source_flux_ratio=1., label='flux ratio = 1',
                         **plot_kwargs)
model.plot_magnification(source_flux_ratio=5., label='flux ratio = 5',
                         **plot_kwargs)
model.plot_magnification(source_flux_ratio=0.2, label='flux ratio = 0.2',
                         **plot_kwargs)
model_1.plot_magnification(label='only source 1', ls='--', **plot_kwargs)
model_2.plot_magnification(label='only source 2', ls='--', **plot_kwargs)

plt.title('Binary-source binary-lens light curves')
plt.legend()

# Make trajectory plots:
plt.figure()
model.plot_trajectory(caustics=True, label='both_sources', **plot_kwargs)
model_2.plot_trajectory(label='source 2', ls='--', lw=4, **plot_kwargs)
model_1.plot_trajectory(label='source 1', ls='--', lw=4, **plot_kwargs)

plt.title('Caustic and trajectories of both sources')
plt.legend()

# Combine it with a dataset (that doesn't show any caustic crossing),
# calculate chi^2, and make a plot of data and model:
file_name = os.path.join(mm.DATA_PATH, "photometry_files", "OB08092",
                         "phot_ob08092_O4.dat")
data = mm.MulensData(file_name=file_name)
event = mm.Event(datasets=data, model=model)
print("chi2 = ", event.get_chi2())
plt.figure()
event.plot_data()
event.plot_model(c='red', **plot_kwargs)
plt.xlim(t_start, t_stop)
plt.show()
