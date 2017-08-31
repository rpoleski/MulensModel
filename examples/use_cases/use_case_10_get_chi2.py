import matplotlib.pyplot as pl
from astropy import units as u 
import os

import MulensModel


data = MulensModel.MulensData(file_name=os.path.join(MulensModel.MODULE_PATH, 
                                                'data', 'phot_ob160023.dat'))

model = MulensModel.Model()
model.set_parameters(t_0=2457518.902, u_0=0.590, t_E=133.34*u.day)

event = MulensModel.Event(datasets=data, model=model) 
# note that data is an instance of Mulens.Data but event.datasets is a list

event.get_chi2()
print(event.chi2)

#############
#Original Use Case
print(event.model.time, event.model.magnification, event.model.flux, 
      event.model.magnitude)

print(event.model.bandpass, event.model.source_flux, event.model.blend_flux) 
# default value of event.model.bandpass is event.data[0].bandpass or 'W149' if no data; for source_flux -> 1.0; for blend_flux -> 0.0

for dataset in event.datasets:
    pl.scatter(dataset.time, dataset.mag)
    pl.scatter(event.model.time_data[dataset], event.model.mag_data[dataset]) 

#############
# Actual Implementation (except that model.flux and model.magnitude
# aren't implemented)

import numpy as np
model_times = np.arange(2457200, 2457800, 100.)

print(
    model_times, 
    event.model.magnification(model_times), 
    event.model.flux(model_times), 
    event.model.magnitude(model_times))

for (i, dataset) in enumerate(event.datasets):
    (f_source, f_blend) = event.model.get_ref_fluxes(data_ref=dataset)

    pl.scatter(dataset.time, dataset.mag, label='DataSet {0}'.format(i))
    pl.scatter(
        dataset.time, 
        event.model.magnitude(dataset.time, f_source=f_source, f_blend=f_blend),
        label='Model for {0}'.format(i))

pl.legend()
pl.show()

#############
#Alternative Use Case (which will require implementation)

print(event.model.time, event.model.magnification, event.model.flux, 
      event.model.magnitude)

print(event.model.bandpass, event.model.source_flux, event.model.blend_flux) 
#Since these weren't set, they should be the default values
#(event.model.bandpass is event.datasets[0].bandpass or 'W149' if no
#data; source_flux and blend_flux are set by the event.datasets[0] or
#source_flux -> 1.0; for blend_flux -> 0.0 if no data

#Were these datasets meant to be scaled to the same flux system?
for (i, dataset) in enumerate(event.datasets):
    pl.scatter(dataset.time, dataset.mag, label='DataSet {0}'.format(i))
    pl.scatter(
        dataset.time, event.model.magnitude(dataset), 
        label='Model for {0}'.format(i)) 

pl.legend()
pl.show()

