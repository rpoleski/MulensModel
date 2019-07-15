import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

sys.path.append("/usr/custom/pyLIMA-1.0.0")
from pyLIMA import event, telescopes, microlmodels
import MulensModel as mm


# Common settings:
t_0 = 2456900.
u_0 = 0.1
t_E = 150.
pi_E_N = 1.
pi_E_E = 2.
ra_deg = 270.
dec_deg = -30.
time = np.linspace(t_0-80., t_0+80., 1000)
d_time = 2450000.

# MulensModel calculations and plot:
coords = SkyCoord(ra_deg, dec_deg, unit=u.deg)
params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E,
          'pi_E_N': pi_E_N, 'pi_E_E': pi_E_E}
mm_model = mm.Model(params, coords=coords)
mm_mag = mm_model.magnification(time)
plt.plot(time-d_time, mm_mag, 'r.', label='MulensModel')

# Read VBBL output and plot it:
vbbl_data = np.loadtxt("fake.out", unpack=True)
plt.plot(vbbl_data[0], vbbl_data[1], 'g-.', label='VBBL')

# This are the changes I have to make to make the results as close as possible:
pi_E_E = -pi_E_E
pi_E_N = -pi_E_N

# pyLIMA calculations and plots:
your_event = event.Event()
your_event.ra = ra_deg
your_event.dec = dec_deg
data_1 = np.array([time, time*0.+15., time*0.+.01]).T
telescope_1 = telescopes.Telescope(light_curve_magnitude=data_1)
your_event.telescopes.append(telescope_1)
your_event.check_event()
model_1 = microlmodels.create_model(
    'PSPL', your_event, parallax=['Annual', t_0])
model_1.define_model_parameters()
pyLIMA_parameters = model_1.compute_pyLIMA_parameters(
    [t_0, u_0, t_E, pi_E_N, pi_E_E])
model = model_1.compute_the_microlensing_model(telescope_1, pyLIMA_parameters)
mag_pyLIMA = model_1.model_magnification(telescope_1, pyLIMA_parameters)
plt.plot(time-d_time, mag_pyLIMA, 'b.', label='pyLIMA')

# Compare pyLIMA and MM:
index = np.argmax(np.abs(mm_mag - mag_pyLIMA))
print("Largest difference is for: ", index, time[index]-d_time)
print("pyLIMA:", mag_pyLIMA[index])
print("MM:", mm_mag[index])

# This is the end:
plt.legend(loc='best')
plt.xlabel('JD-2450000')
plt.ylabel('magnification')
plt.show()
