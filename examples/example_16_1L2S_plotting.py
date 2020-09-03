"""
This example checks plotting functions for binary sources. It is
derived from example_11_binary_source.py

New code starts at Line 38.
"""
import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np

# First, prepare the data. There is nothing very exciting in this part,
# so you may skip it.
t_0_1 = 6100.
u_0_1 = 0.2
t_0_2 = 6140.
u_0_2 = 0.01
t_E = 25.
assumed_flux_1 = 100.
assumed_flux_2 = 5.
assumed_flux_blend = 10.
n_a = 1000
n_b = 600
time_a = np.linspace(6000., 6300., n_a)
time_b = np.linspace(6139., 6141., n_b)
time = np.sort(np.concatenate((time_a, time_b)))
model_1 = mm.Model({'t_0': t_0_1, 'u_0': u_0_1, 't_E': t_E})
A_1 = model_1.magnification(time)
model_2 = mm.Model({'t_0': t_0_2, 'u_0': u_0_2, 't_E': t_E})
A_2 = model_2.magnification(time)
flux = A_1 * assumed_flux_1 + A_2 * assumed_flux_2 + assumed_flux_blend
flux_err = 6. + 0. * time
flux += flux_err * np.random.normal(size=n_a+n_b)
my_dataset = mm.MulensData([time, flux, flux_err], phot_fmt='flux')

# Model
params = {'t_0_1': 6101., 'u_0_1': 0.19, 't_0_2': 6140.123, 'u_0_2': 0.04,
          't_E': 20.}
my_model = mm.Model(params)
my_event = mm.Event(datasets=my_dataset, model=my_model)

# NEW CODE STARTS HERE
# Plot the model from Event()
plt.figure()
plt.title('Event() Model')
my_event.plot_model()

(source_flux, blend_flux) = my_event.get_ref_fluxes()
# Plot just the model
plt.figure(figsize=(6, 8))
# Plot the model in "effective" magnification
plt.subplot(3, 1, 1)
plt.title('Model Magnification')
my_model.plot_magnification(q_flux=assumed_flux_2/assumed_flux_1)
# Plot the model in magnitudes
# specifying f_source
plt.subplot(3, 1, 2)
plt.title('Model Lightcurve with Specified Fluxes')
my_model.plot_lc(f_source=source_flux, f_blend=blend_flux)
# specifying q_flux (using assumed values, so should be different from prev,
# which uses fitted fluxes)
plt.subplot(3, 1, 3)
plt.title('Model Lightcurve with q_flux')
my_model.plot_lc(
    f_source=assumed_flux_1, f_blend=assumed_flux_blend,
    q_flux=assumed_flux_2/assumed_flux_1)
plt.tight_layout()

plt.show()