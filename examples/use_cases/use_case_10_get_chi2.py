import matplotlib.pyplot as pl
from astropy import units as u 

import MulensModel


data=MulensModel.MulensData(file_name='some_file.dat')

model = MulensModel.Model()
model.parameters(t_0=7600., u_0=0.01, t_E=34.*u.day)

e = MulensModel.Event(datasets=data, model = model) # note that data is an instance of Mulens.Data but e.datasets is a list

e.get_chi2()

print(e.model.time, e.model.A, e.model.flux, e.model.mag)

print(e.chi2)

print(e.model.bandpass, e.model.source_flux, e.model.blend_flux) # default value of e.model.bandpass is e.data[0].bandpass or 'W149' if no data; for source_flux -> 1.0; for blend_flux -> 0.0

for dataset in e.datasets:
    pl.plot(dataset.time, dataset.mag)
    pl.plot(e.model.time_data[dataset], e.model.mag_data[dataset]) 


