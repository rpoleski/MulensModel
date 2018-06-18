import os
import sys
import numpy as np

sys.path.append("/usr/custom/pyLIMA-1.0.0") 
from pyLIMA import event, telescopes, microlmodels
from pyLIMA_functions import chi2_telescope


file_1 = "test_10000.txt"
data_1 = np.loadtxt(file_1)
telescope_1 = telescopes.Telescope(light_curve_magnitude=data_1)

your_event = event.Event()
your_event.telescopes.append(telescope_1)
model_1 = microlmodels.create_model('PSPL', your_event)
model_1.define_model_parameters()

t_0 = 2456900.
u_0 = 0.01
t_E = 20.

parameters_list = [t_0, u_0, t_E]
