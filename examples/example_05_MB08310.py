"""
Plots data and model for MB08310 with residuals (no fitting). 

From Janczak et al. 2010, ApJ 711, 731
"""
import glob
import os
import matplotlib.pyplot as pl
from matplotlib import gridspec

import MulensModel
from MulensModel.mulensdata import MulensData
from MulensModel.event import Event
from MulensModel.model import Model

#Read in MB08310 data files (see data/MB08310) as MulensData objects
MODULE_PATH = "/".join(MulensModel.__file__.split("/source")[:-1])
print('MODULE_PATH: {0}'.format(MODULE_PATH))

#Grabbing all data files in the MB08310 folder
files = glob.glob(MODULE_PATH+"/data/MB08310"+"/*.tbl")
datasets = []
labels = []
for file_ in sorted(files):
    data = MulensData(file_name=file_, comments=["\\","|"])	
    datasets.append(data)
    labels.append(os.path.basename(file_))

#Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
rho = t_star / t_E
plens_model = Model(t_0=t_0, u_0=u_0, t_E=t_E, rho=rho)
plens_model.set_magnification_methods([t_0-.05, 'finite_source_uniform_Gould94', t_0+.05])

#Combine the data and model into an event
ev = Event(datasets=datasets, model=plens_model)
ev.data_ref = 6

gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])

#Plot the data and model
pl.figure()
pl.subplot(gs[0])
ev.plot_model(subtract_2450000=True)
ev.plot_data(subtract_2450000=True)
pl.title('Data and Fitted Model (Default)')
#Plot the residuals 
pl.subplot(gs[1])
ev.plot_residuals(subtract_2450000=True)

#Plot the data and model (customized)
pl.figure()
pl.subplot(gs[0])
t_start= t_0 - 3.
t_stop = t_0 + 1.
ev.plot_model(color='black', t_start=t_start, t_stop=t_stop, subtract_2450000=True)
ev.plot_data(
    label_list=labels, marker='o', markersize=5,  
    color_list=['black', 'red', 'yellow', 'green', 'cyan', 'blue', 'purple'],
    subtract_2450000=True)
pl.ylim(17.5, 12.5)
pl.xlim(t_start-2450000., t_stop-2450000.)
pl.legend(loc='upper left')
pl.title('Data and Fitted Model (Custom)')

#Plot the residuals 
pl.subplot(gs[1])
ev.plot_residuals(subtract_2450000=True)
pl.xlim(t_start-2450000., t_stop-2450000.)

pl.show()
