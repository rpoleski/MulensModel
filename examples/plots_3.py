"""
Figure 3.

How to create and plot a model and then add some data and fit for the
source and blend fluxes.

This example shows OGLE-2003-BLG-235/MOA-2003-BLG-53, the first
microlensing planet. See [Bond et
al. 2004](http://adsabs.harvard.edu/abs/2004ApJ...606L.155B). The data
were downloaded from the [NASA Exoplanet
Archive](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=OGLE-2003-BLG-235L+b&type=CONFIRMED_PLANET).
"""

# Import basic packages
import MulensModel
import matplotlib.pyplot as pl # MulensModel uses matplotlib for plotting.

# Import data
OGLE_data = MulensModel.MulensData(
    file_name='../data/OB03235/OB03235_OGLE.tbl.txt', comments=['\\','|'])
MOA_data = MulensModel.MulensData(
    file_name='../data/OB03235/OB03235_MOA.tbl.txt', phot_fmt='flux', 
    comments=['\\','|'])

# Define a model with 2-bodies:
my_1S2L_model = MulensModel.Model(
    {'t_0': 2452848.06, 'u_0': 0.133, 't_E': 61.5, 'rho': 0.00096, 
     'q': 0.0039, 's': 1.120, 'alpha': 223.8})
# Since rho is set, define a time range and method to apply 
# finite source effects:
my_1S2L_model.set_magnification_methods(
    [2452833., 'VBBL', 2452845.])

# Combine the data and model into an Event:
my_event = MulensModel.Event(
    datasets=[MOA_data, OGLE_data], model=my_1S2L_model)

# Plot the result:
my_event.plot_model(
    t_range=[2452810,2452890], subtract_2450000=True, color='black',
    data_ref=1)
my_event.plot_data(
    subtract_2450000=True, data_ref=1, label_list=['MOA', 'OGLE'], 
    color_list=['cyan', 'orange'], s=5)
# MulensModel automatically fits for the source and blend flux for the  
# given model.

# Customize the output
pl.legend(loc='best')
pl.title('OGLE-2003-BLG-235/MOA-2003-BLG-53')
pl.ylim(19., 16.5)
pl.xlim(2810,2890)
pl.savefig('figure_3.png')
