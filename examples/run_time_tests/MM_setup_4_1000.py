import os
import sys
sys.path.append("/home/antares/poleski/WFIRST/MulensModel/source")
import MulensModel as mm


file_1 = "test_1000_piE.txt"
data_1 = mm.MulensData(file_name=file_1)
t_0 = 2456900.
u_0 = 0.1
t_E = 50.
pi_E_N = 0.6
pi_E_E = 0.8
param = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E_N': pi_E_N, 'pi_E_E': pi_E_E, 't_0_par': t_0}
model = mm.Model(param, coords="18:00:00.00 -30:00:00.0")
event = mm.Event(datasets=[data_1], model=model)
event.sum_function = 'numpy.sum'

event.get_chi2()
