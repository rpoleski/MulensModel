import os
import sys
sys.path.append("/home/antares/poleski/WFIRST/MulensModel/source")
import MulensModel as mm


file_1 = "/home/antares/poleski/WFIRST/MulensModel/data/photometry_files/phot_ob08092_O4.dat"
data_1 = mm.MulensData(file_name=file_1)
data_1.input_fmt = 'flux'
t_0 = 5379.57091
u_0 = 0.52298
t_E = 17.94002
model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
event = mm.Event(datasets=[data_1], model=model)
