import os
import matplotlib.pyplot as pl
from astropy import units as u

import MulensModel


#Read in some data
data = []
data.append(MulensModel.MulensData(file_name=os.path.join(
                MulensModel.MODULE_PATH, 'data', 'ob151100_OGLE_ref_v1.dat'))

pl.figure()
pl.errorbar(data[0].time, data[0].mag, yerr=data[0].err_mag, fmt='o')
pl.title('Raw Data (MulensData)')
pl.gca().invert_yaxis()

#Generate a model
model = MulensModel.Model()
model.set_parameters(t_0=2457181.9, u_0=0.088, t_E=20.291*u.day)

pl.figure()
model.plot_lc(f_source=1.0, f_blend=0.0)
pl.title('Base Model')

#Combine Model and Data and plot
event = MulensModel.Event(datasets=data, model=model)

#Plot the data
pl.figure()
pl.subplot(2, 1, 1)
event.plot_data()
event.plot_model()
pl.title('Data and Model')

pl.subplot(2, 1, 2)
event.plot_residuals()

pl.show()

