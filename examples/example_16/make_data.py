"""
Prepare the ob008092 photometry with added periodic signal
"""
import numpy as np

import MulensModel as mm


file_in = "data/OB08092/phot_ob08092_O4.dat"
file_out = "data/OB08092/phot_ob08092_O4_periodic.dat"
t_0 = 5500.
period = 12.34567
amplitude = 100


data_in = mm.MulensData(file_name=file_in)

flux = data_in.flux + amplitude * np.sin(2*np.pi*(data_in.time - t_0)/period)

mag = mm.Utils.get_mag_from_flux(flux)

out = ""
for (t, m, e) in zip(data_in.time, mag, data_in.err_mag):
    out += "{:.5f} {:.4f} {:.4f}\n".format(t, m, e)

with open(file_out, "w") as file_out:
    file_out.write(out)
