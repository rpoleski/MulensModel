"""
This example checks plotting functions for binary sources. It is
derived from example_11_binary_source.py
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
model_1 = mm.Model({'t_0': t_0_1, 'u_0': u_0_1, 't_E': t_E})
model_2 = mm.Model({'t_0': t_0_2, 'u_0': u_0_2, 't_E': t_E})


def generate_time_vector(n_a, n_b):
    """Generate sorted array simulating survey + follow-up time vector."""
    time_a = np.linspace(6000., 6300., n_a)
    time_b = np.linspace(6139., 6141., n_b)
    time = np.sort(np.concatenate((time_a, time_b)))
    return time


def generate_dataset(time, flux_1, flux_2, blend_flux, flux_err,
                     model_1, model_2):
    """Generate simulated dataset assuming binary source model."""
    A_1 = model_1.get_magnification(time)
    A_2 = model_2.get_magnification(time)
    flux = A_1 * flux_1 + A_2 * flux_2 + blend_flux
    err_flux = flux_err + 0. * time
    flux += flux_err * np.random.normal(size=len(time))
    my_dataset = mm.MulensData([time, flux, err_flux], phot_fmt='flux')
    return my_dataset


assumed_flux_1 = 100.
assumed_flux_2 = 5.
assumed_flux_blend = 10.
n_a = 1000
n_b = 600
flux_err = 6.

time = generate_time_vector(n_a, n_b)
my_dataset = generate_dataset(
    time, assumed_flux_1, assumed_flux_2, assumed_flux_blend, flux_err,
    model_1, model_2)

time_2 = generate_time_vector(int(n_a / 5), int(n_b / 5))
my_dataset_2 = generate_dataset(
    time_2, assumed_flux_1/2., assumed_flux_2/2., assumed_flux_blend/2.,
    2.*flux_err, model_1, model_2)

# Model
params = {'t_0_1': t_0_1, 'u_0_1': u_0_1, 't_0_2': t_0_2, 'u_0_2': u_0_2,
          't_E': t_E}
my_model = mm.Model(params)
my_event = mm.Event(datasets=[my_dataset, my_dataset_2], model=my_model)

# Plot just the data
plt.figure()
plt.title('Raw Data')
my_dataset.plot(phot_fmt='mag')
my_dataset_2.plot(phot_fmt='mag')
(source_flux, blend_flux) = my_event.get_ref_fluxes()

# Plot just the model:
# Plot the model in "effective" magnification
plt.figure(figsize=(6, 8))
plt.subplot(3, 1, 1)
plt.title('Model Magnification')
my_model.plot_magnification(source_flux_ratio=assumed_flux_2 / assumed_flux_1)
# Plot the model in magnitudes
# - specifying f_source
plt.subplot(3, 1, 2)
plt.title('Model Lightcurve with Specified Fluxes')
my_model.plot_lc(source_flux=source_flux, blend_flux=blend_flux)
# - specifying q_flux (using assumed values, so should be different from prev,
#   which uses fitted fluxes)
plt.subplot(3, 1, 3)
plt.title('Model Lightcurve with q_flux')
my_model.plot_lc(
    source_flux=assumed_flux_1, blend_flux=assumed_flux_blend,
    source_flux_ratio=assumed_flux_2 / assumed_flux_1)
plt.tight_layout()

# Plot the model and data from Event()
plt.figure()
plt.title('Event() Model + Data')
my_event.plot_model(zorder=10, color='black')
my_event.plot_data()

# Plot the source trajectories
plt.figure()
plt.title('Model Source Trajectories')
my_model.plot_trajectory()

plt.show()
