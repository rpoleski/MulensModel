"""
Plots data and model for MB08310 with residuals (no fitting).

From `Janczak et al. 2010, ApJ 711, 731
<https://ui.adsabs.harvard.edu/abs/2010ApJ...711..731J/abstract>`_.

"""
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

import MulensModel as mm


# Read in MB08310 data files (see data/MB08310) as MulensData objects.
# Grabbing all data files in the MB08310 folder
files = glob.glob(os.path.join(mm.DATA_PATH, "photometry_files",
                               "MB08310", "*.tbl"))

datasets_default = []
for file_ in sorted(files):
    data = mm.MulensData(file_name=file_, comments=["\\", "|"])
    datasets_default.append(data)

# Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
plens_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star})
method = 'finite_source_uniform_Gould94'
plens_model.set_magnification_methods([t_0-2.*t_star, method, t_0+2.*t_star])

# Combine the data and model into an event
event_default = mm.Event(datasets=datasets_default, model=plens_model)
event_default.data_ref = 6

# F1 & F2: (data and model) and trajectory using minimal commands.
event_default.plot(
    subtract_2450000=True, trajectory=True, title='Minimal Effort Plot')

# F3: Plot the data and model
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
plt.figure()
ax11 = plt.subplot(gs[0])
event_default.plot_model(subtract_2450000=True)
event_default.plot_data(subtract_2450000=True)
plt.title('Data and Fitted Model (Default)')
# Plot the residuals
plt.subplot(gs[1], sharex=ax11)
event_default.plot_residuals(subtract_2450000=True)

# F4: Plot the trajectory
plt.figure()
plt.title('Trajectory w/Data (Default)')
event_default.plot_trajectory()
event_default.plot_source_for_datasets()

# -----------------
# F5: Plot the data and model (customized)
datasets_custom = []
color_list = ['black', 'red', 'yellow', 'green', 'cyan', 'blue', 'purple']
for (color, file_) in zip(color_list, sorted(files)):
    label = os.path.basename(file_).split('_', maxsplit=2)[0]
    if label == 'CTIO':
        label += ' ' + os.path.basename(file_).split('_', maxsplit=2)[1]

    data = mm.MulensData(
        file_name=file_, comments=["\\", "|"],
        plot_properties={
            'color': color,
            'label': label})

    datasets_custom.append(data)

event_custom = mm.Event(datasets=datasets_custom, model=plens_model)

t_start = t_0 - 3.
t_stop = t_0 + 1.
n_star = 2.
t_star_start = t_0 - n_star * t_star
t_star_stop = t_0 + n_star * t_star
times = np.arange(t_start, t_star_start, 0.01)
times = np.concatenate((times, np.arange(t_star_start, t_star_stop, 0.0001)))
times = np.concatenate((times, np.arange(t_star_stop, t_stop, 0.01)))

plt.figure()
ax31 = plt.subplot(gs[0])
event_custom.plot_model(
    times=times, color='black', subtract_2450000=True)
event_custom.plot_data(marker='s', markersize=3, subtract_2450000=True)
plt.ylim(17.5, 12.5)
plt.xlim(t_start-2450000., t_stop-2450000.)
plt.legend(loc='upper left')
plt.title('Data and Fitted Model (Custom)')

# Plot the residuals
plt.subplot(gs[1], sharex=ax31)
event_custom.plot_residuals(marker='s', markersize=3, subtract_2450000=True)
plt.xlim(t_start-2450000., t_stop-2450000.)

# F6: Plot the trajectory
plt.figure()
plt.title('Trajectory w/Data (Custom)')
plt.gca().set_aspect('equal')
event_custom.plot_trajectory()
event_custom.plot_source_for_datasets()
traj_range = (-0.05, 0.05)
plt.xlim(traj_range)
plt.ylim(traj_range)

plt.show()
