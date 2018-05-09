import MulensModel as mm
import matplotlib.pyplot as pl
import numpy as np

raise NotImplementedError('This use case has not been implemented.')

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
times = np.union1d(
    np.arange(t_0_1 - n_t_E * t_E, t_0_1 + n_t_E * t_E, 2. * n_t_E * t_E / n_pts),
    np.arange(t_0_2 - u_0_2 * t_E, t_0_2 + u_0_2 * t_E, 2. * u_0_2 * t_E / 100.))

# Distinct Models for the two sources.
model_1 = mm.Model(dict(t_0=t_0_1, u_0=u_0_1, t_E=t_E))
model_2 = mm.Model(dict(t_0=t_0_2, u_0=u_0_2, t_E=t_E))

combined_model_lc = mm.Utils.get_mag_from_flux(
    f_source_1 * model_1.magnification(times) + f_source_2 * model_2.magnification(times) + f_blend)

# Binary Source Model version
binary_source_model = mm.Model(dict(
   t_0_1=t_0_1, u_0_1=u_0_1, t_0_2=t_0_2, u_0_2=u_0_2, t_E=t_E, flux_ratio=f_source_2/f_source_1))

# Plots
pl.figure()
pl.title('2 Point Source Models')
pl.plot(times, combined_model_lc, color='blue')
pl.gca().invert_yaxis()

pl.figure()
pl.title('1 Binary Source Model')
binary_source_model.plot_lc(times, f_source=f_source_1, color='magenta')

pl.show()