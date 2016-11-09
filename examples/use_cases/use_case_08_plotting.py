import matplotlib.pyplot as pl
from astropy import units as u

import MulensModel

data = []
data.append(MulensModel.MulensData(file_name = 'some_file.dat'))
pl.plot(data[0].time('JD'), data[0].mag)
pl.plot(data[0].hjd, data[0].mag)

model = MulensModel.Model()
model.parameters(t_0=7600., u_0=0.01, t_E=34.*u.day)
pl.plot(model.time, model.magnitude) # Symmetry with data plotting?

event = MulensModel.Event(datasets=data, model=model)

### jcy - I find this implementation dissatisfying
pl.subplot(2, 1, 1)
pl.plot(event.fit.time[0], event.fit.magnitude_data[0])
pl.plot(event.fit.time[0], event.fit.magnitude_model[0])

pl.subplot(2, 1, 2)
pl.plot(event.fit.time[0], event.fit.res[0]) #res_data?

