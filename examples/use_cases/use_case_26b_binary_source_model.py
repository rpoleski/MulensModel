import matplotlib.pyplot as plt
import numpy as np

import MulensModel as mm

# Model parameters
t_0_1 = 2458000.
u_0_1 = 0.1
t_0_2 = 2458020.
u_0_2 = 0.01
t_E = 30.

# Source fluxes
f_source_1 = 1.
f_source_2 = 0.01
f_blend = 0.

# Times to generate the model
n_t_E = 3.
n_pts = 1000
dt_1 = n_t_E * t_E
dt_2 = u_0_2 * t_E
times = np.union1d(
    np.arange(t_0_1 - dt_1, t_0_1 + dt_1, 2. * dt_1 / n_pts),
    np.arange(t_0_2 - dt_2, t_0_2 + dt_2, 2. * dt_2 / 100.))

# Distinct Models for the two sources.
model_1 = mm.Model(dict(t_0=t_0_1, u_0=u_0_1, t_E=t_E))
model_2 = mm.Model(dict(t_0=t_0_2, u_0=u_0_2, t_E=t_E))

mag_1 = model_1.get_magnification(times)
mag_2 = model_2.get_magnification(times)
combined_model_lc = mm.Utils.get_mag_from_flux(
    f_source_1 * mag_1 + f_source_2 * mag_2 + f_blend)

# Binary Source Model version
binary_source_model = mm.Model(dict(
    t_0_1=t_0_1, u_0_1=u_0_1, t_0_2=t_0_2, u_0_2=u_0_2, t_E=t_E))

# Plots
plt.figure()
plt.title('Sum of Two Point Source Models')
plt.plot(times, combined_model_lc, color='blue')
plt.gca().invert_yaxis()

plt.figure()
plt.title('Binary Source Model')
binary_source_model.plot_lc(
    times, source_flux=f_source_1, source_flux_ratio=f_source_2/f_source_1,
    blend_flux=f_blend, color='magenta')

plt.show()
