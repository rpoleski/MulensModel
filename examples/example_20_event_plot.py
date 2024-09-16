"""
Show  (and test) options for event.plot().
"""
import os.path
import matplotlib.pyplot as plt
import glob

import MulensModel as mm

# PSPL Event (Derived from Example 05)
files = glob.glob(os.path.join(mm.DATA_PATH, "photometry_files",
                               "MB08310", "*.tbl"))

# Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
pspl_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star})
method = 'finite_source_uniform_Gould94'
pspl_model.set_magnification_methods([t_0-2.*t_star, method, t_0+2.*t_star])

datasets_custom = []
color_list = ['black', 'red', 'yellow', 'green', 'cyan', 'blue', 'purple']
for (color, file_) in zip(color_list, sorted(files)):
    data = mm.MulensData(
        file_name=file_, comments=["\\", "|"],
        plot_properties={
            'color': color,
            'label': os.path.basename(file_).split('_', maxsplit=2)[0]})
    datasets_custom.append(data)

pspl_event = mm.Event(datasets=datasets_custom, model=pspl_model)
pspl_t_start = t_0 - 3.
pspl_t_stop = t_0 + 1.

# Default Behavior
# Fig 1: Light curve with data, model, and residuals
pspl_event.plot()

# With options
# Fig 2: Light curve with data and model, no residuals, no errorbars.
# Fig 3: Source trajectory with points plotted for each observation
pspl_event.plot(
    t_range=[pspl_t_start, pspl_t_stop], residuals=False, show_errorbars=False,
    trajectory=True, title='MB08310', subtract_2450000=False)

# Planetary Lens Event (Derived from Example 16)
dir_1 = os.path.join(mm.DATA_PATH, "photometry_files", 'OB03235')
file_1 = os.path.join(dir_1, 'OB03235_OGLE.tbl.txt')
ogle_data = mm.MulensData(
    file_name=file_1, bandpass='I', comments=["\\", "|"],
    plot_properties={'zorder': 10., 'color': 'red', 'label': "OGLE I-band"})
file_2 = os.path.join(dir_1, 'OB03235_MOA.tbl.txt')
moa_data = mm.MulensData(
        file_name=file_2, phot_fmt='flux', comments=["\\", "|"],
        plot_properties={'show_errorbars': False})
planet_datasets = [ogle_data, moa_data]
planet_model = mm.Model(
    {'t_0': 2452848.06, 'u_0': 0.133, 't_E': 61.5, 'rho': 0.00096,
     'q': 0.0039, 's': 1.120, 'alpha': 43.8})
planet_model.set_magnification_methods([2452833., 'VBBL', 2452845.])
planet_event = mm.Event(
    datasets=planet_datasets, model=planet_model)

# Default Behavior
# Fig 4: Light curve with data, model, and residuals. No MOA errorbars.
# Fig 5: Source trajectory with points plotted for each observation.
planet_event.plot()

# With options
# Fig 6: Light curve with data, model, and residuals. Errorbars for all data
# points. No trajectory plot.
planet_event.plot(
    t_range=[2452810, 2452890], residuals=True, show_errorbars=True,
    legend=False, trajectory=False, title='OB03235',
    subtract_2450000=True, subtract_2460000=False, data_ref=0)

plt.show()
