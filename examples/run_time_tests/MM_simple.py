import MulensModel as mm
import os


data_path = os.path.abspath(__file__)
for i in range(3):
    data_path = os.path.dirname(data_path)
data_path = os.path.join(data_path, 'data') 

file_1 = os.path.join(data_path, "photometry_files", "phot_ob08092_O4.dat")
data_1 = mm.MulensData(file_name=file_1)

t_0 = 5379.57091
u_0 = 0.52298
t_E = 17.94002

model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
event = mm.Event(datasets=[data_1], model=model)

print(event.get_chi2())

