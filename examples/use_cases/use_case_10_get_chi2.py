import matplotlib.pyplot as pl
from astropy import units as u 

import MulensModel


data=MulensModel.MulensData(file_name='some_file.dat')

model = MulensModel.Model()
model.parameters(t_0=7600., u_0=0.01, t_E=34.*u.day)

event = MulensModel.Event(datasets=data, model=model) 
# note that data is an instance of Mulens.Data but event.datasets is a list

event.get_chi2()

print(event.model.time, event.model.magnification, event.model.flux, 
      event.model.magnitude)

print(event.chi2)

print(event.model.bandpass, event.model.source_flux, event.model.blend_flux) 
# default value of event.model.bandpass is event.data[0].bandpass or 'W149' if no data; for source_flux -> 1.0; for blend_flux -> 0.0

for dataset in event.datasets:
    pl.plot(dataset.time, dataset.mag)
    pl.plot(event.model.time_data[dataset], event.model.mag_data[dataset]) 


