import os
import sys
import numpy as np

sys.path.append("/usr/custom/pyLIMA-1.0.0") 
from pyLIMA import event, telescopes, microlmodels
from pyLIMA_functions import chi2_telescope

t_0 = 2456900.
u_0 = 0.1
t_E = 50.
pi_E_N = 0.6
pi_E_E = 0.8

pi_E_N = -pi_E_N
pi_E_E = -pi_E_E

file_1 = "test_1000_piE.txt"
data_1 = np.loadtxt(file_1)
telescope_1 = telescopes.Telescope(light_curve_magnitude=data_1)

your_event = event.Event()
your_event.ra = 270.
your_event.dec = -30.
your_event.telescopes.append(telescope_1)
model_1 = microlmodels.create_model('PSPL', your_event, parallax=['Annual', t_0])
model_1.define_model_parameters()

parameters_list = [t_0, u_0, t_E, pi_E_N, pi_E_E]

chi2_telescope(your_event, model_1, parameters_list)
