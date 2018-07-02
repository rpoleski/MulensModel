import os
import sys
sys.path.append("/home/antares/poleski/WFIRST/MulensModel/source")
import MulensModel as mm


file_1 = "test_100.txt"
data_1 = mm.MulensData(file_name=file_1)
t_0 = 2456900.
u_0 = 0.01
t_E = 20.
model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
event = mm.Event(datasets=[data_1], model=model)
event.sum_function = 'numpy.sum'

