"""
Plot a binary lens, triple source model with trajectory and caustics including orbital motion and parallax.

Also show the __repr__ string for that model.

Uses OB151459 as a template.
"""
import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np

# 2L2S event
subhead_2L2S = 'Table 5: 2L2S Model - close '

coords = '18:00:50.40 −28:40:15.7'

t_E = 4.613
s = 0.857
q = 2.374e-3
alpha = -np.rad2deg(0.890)
params_static = {
    't_E': t_E, 's': s, 'q': q, 'alpha': alpha,
    't_0_1': 7199.941, 'u_0_1': 0.069, 'rho_1': 0.688e-3,
    't_0_2': 7200.057, 'u_0_2': 0.039, 'rho_2': 0.542e-3,
    't_0_3': 7197., 'u_0_3': -0.28, 'rho_3': 2.e-3}

params = {**params_static,
    'pi_E_E': 100.0, 'pi_E_N': 1.0,
    'ds_dt': 3.0, 'dalpha_dt': 360., 't_0_kep': 7200.2}

q_F_2 = 0.181
source_flux_1 = 0.061
source_flux_2 = source_flux_1 * q_F_2
source_flux_3 = 3. * source_flux_1
blend_flux = 0.095
source_fluxes = [source_flux_1, source_flux_2, source_flux_3]

mag_methods_12 = [7200.1, 'VBBL', 7200.3]
mag_methods_3 = [7195., 'VBBL', 7197.5]

model = mm.Model(params, coords=coords)
model.set_magnification_methods(mag_methods_12, 1)
model.set_magnification_methods(mag_methods_12, 2)
model.set_magnification_methods(mag_methods_3, 3)
print(model)

model_static = mm.Model(params_static, coords=coords)
model_static.set_magnification_methods(mag_methods_12, 1)
model_static.set_magnification_methods(mag_methods_12, 2)
model_static.set_magnification_methods(mag_methods_3, 3)

# Make the plots
times_perturb_1 = np.arange(7200.1, 7200.3, 0.001)
times_perturb_2 = np.arange(7200.195, 7200.21, 0.0001)
times_lc = np.arange(7190., 7210., 0.01)
times = np.hstack((times_lc, times_perturb_1, times_perturb_2))
times.sort()

plt.figure()
model.plot_trajectory(caustics=True)
model.plot_caustics(epoch=7196., color='blue')
plt.gca().set_aspect('equal')
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.minorticks_on()

plt.figure()
model_static.plot_lc(
    times, label='static', color='gray', linestyle=':', source_flux=source_fluxes, blend_flux=blend_flux)
model.plot_lc(times, label='full', source_flux=source_fluxes, blend_flux=blend_flux)
plt.legend()
plt.minorticks_on()
plt.tight_layout()

plt.show()