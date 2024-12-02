"""
Creates last figure of Poleski and Yee (2019)
"Modeling microlensing events with MulensModel"
Astronomy and Computing 26, 35
https://ui.adsabs.harvard.edu/abs/2019A&C....26...35P/abstract
https://arxiv.org/abs/1803.01003

This example shows OGLE-2003-BLG-235/MOA-2003-BLG-53,
the first microlensing planet.

"""
from matplotlib import pyplot
import os

from MulensModel import MulensData, Model, Event, MODULE_PATH


# Import data
data_dir = os.path.join(MODULE_PATH, 'data', 'photometry_files', 'OB03235')

OGLE_data = MulensData(
    file_name=os.path.join(data_dir, 'OB03235_OGLE.tbl.txt'),
    comments=['\\', '|'])
MOA_data = MulensData(
    file_name=os.path.join(data_dir, 'OB03235_MOA.tbl.txt'),
    comments=['\\', '|'], phot_fmt='flux')

# Define a model with a 2-body lens (these parameters slightly differ
# from Bond et al. 2004):
model_1S2L = Model({'t_0': 2452848.06, 'u_0': 0.1317, 't_E': 61.5,
                    'rho': 0.00096, 'q': 0.0039, 's': 1.120, 'alpha': 43.72})

# Since rho is set, define a time range and method to apply finite
# source effects:
model_1S2L.set_magnification_methods([2452833., 'VBBL', 2452845.])

# Combine the data and model into an Event:
my_event = Event(datasets=[OGLE_data, MOA_data], model=model_1S2L)

# Make the plot:
t_range = [2452800., 2452875.]
pyplot.axes([0.09, 0.08, 0.9, 0.9])
my_event.plot_data(
    subtract_2450000=True, label_list=['OGLE', 'MOA'],
    color_list=['red', 'blue'], zorder_list=[2, 1], s=6)
my_event.plot_model(
    subtract_2450000=True, t_range=t_range, n_epochs=4000, color='black')

pyplot.legend(loc='best')
pyplot.xlim(t_range[0]-2450000., t_range[1]-2450000.)
pyplot.ylim(19.0, 16.7)

pyplot.axes([0.17, 0.7, 0.3, 0.2])  # Figure inset stars here.
model_1S2L.plot_trajectory(caustics=True)
pyplot.xlim(-0.1, 0.45)
pyplot.ylim(-0.14, 0.14)
pyplot.savefig('figure_4.png')
