"""
Plots data and model for MB08310 with residuals (no fitting).

From `Janczak et al. 2010, ApJ 711, 731
<https://ui.adsabs.harvard.edu/abs/2010ApJ...711..731J/abstract>`_.

"""
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec

from MulensModel import Event, Model, MulensData, DATA_PATH

# Read in MB08310 data files (see data/MB08310) as MulensData objects.
# Grabbing all data files in the MB08310 folder
files = glob.glob(os.path.join(DATA_PATH, "photometry_files",
                               "MB08310", "*.tbl"))

datasets_default = []
for file_ in sorted(files):
    data = MulensData(file_name=file_, comments=["\\", "|"])
    datasets_default.append(data)

# Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
plens_model = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star})
method = 'finite_source_uniform_Gould94'
plens_model.set_magnification_methods([t_0-.05, method, t_0+.05])

# Combine the data and model into an event
event_default = Event(datasets=datasets_default, model=plens_model)
event_default.data_ref = 6

gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])

# Plot the data and model
plt.figure()
plt.subplot(gs[0])
event_default.plot_model(subtract_2450000=True)
event_default.plot_data(subtract_2450000=True)
plt.title('Data and Fitted Model (Default)')
# Plot the residuals
plt.subplot(gs[1])
event_default.plot_residuals(subtract_2450000=True)

# -----------------
# Plot the data and model (customized)
datasets_custom = []
color_list = ['black', 'red', 'yellow', 'green', 'cyan', 'blue', 'purple']
for (i, file_) in enumerate(sorted(files)):
    data = MulensData(
        file_name=file_, comments=["\\", "|"],
        plot_properties={
            'color': color_list[i],
            'label': os.path.basename(file_).split('_', maxsplit=2)[0]})
    datasets_custom.append(data)

event_custom = Event(datasets=datasets_custom, model=plens_model)

plt.figure()
plt.subplot(gs[0])
t_start = t_0 - 3.
t_stop = t_0 + 1.
event_custom.plot_model(
    color='black', t_start=t_start, t_stop=t_stop, subtract_2450000=True)
event_custom.plot_data(marker='s', markersize=3, subtract_2450000=True)
plt.ylim(17.5, 12.5)
plt.xlim(t_start-2450000., t_stop-2450000.)
plt.legend(loc='upper left')
plt.title('Data and Fitted Model (Custom)')

# Plot the residuals
plt.subplot(gs[1])
event_custom.plot_residuals(marker='s', markersize=3, subtract_2450000=True)
plt.xlim(t_start-2450000., t_stop-2450000.)

plt.show()
