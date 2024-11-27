"""
Plot model, data, and model together with data and residuals
"""
import os
import matplotlib.pyplot as plt

import MulensModel as mm


# Read in some data
data = []
file_name = os.path.join(
    mm.DATA_PATH, 'photometry_files', 'OB140939', 'ob140939_OGLE.dat')
data.append(mm.MulensData(file_name=file_name))

plt.figure()
plt.errorbar(data[0].time-2450000., data[0].mag, yerr=data[0].err_mag, fmt='o')
plt.title('Raw Data (MulensData)')
plt.gca().invert_yaxis()

# Generate a model
model = mm.Model({'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87})

plt.figure()
model.plot_lc(source_flux=1.0, blend_flux=0.0, subtract_2450000=True)
plt.title('Base Model')

# Combine Model and Data and plot
event = mm.Event(datasets=data, model=model)

# Plot the data
plt.figure()
plt.subplot(2, 1, 1)
event.plot_data(subtract_2450000=True)
event.plot_model(subtract_2450000=True)
plt.title('Data and Model')

plt.subplot(2, 1, 2)
event.plot_residuals(subtract_2450000=True)

plt.show()
