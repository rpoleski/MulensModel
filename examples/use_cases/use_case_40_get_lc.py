"""
Use case showing API for Event.get_scaled_lc() and Event.plot_lc()
"""
import matplotlib.pyplot as plt
import os

import MulensModel as mm


raise NotImplementedError(
    "Event.get_scaled_lc(), Event.plot_lc() and satellite_skycoord keyword in "
    "Model.get_scaled_lc()/plot_lc() are not implemented yet")

# Import data
dir_ = os.path.join(mm.DATA_PATH, "photometry_files", "OB140939")
file_name_1 = os.path.join(dir_, "ob140939_OGLE.dat")
file_name_2 = os.path.join(dir_, "ob140939_Spitzer.dat")
ephemeris_2 = os.path.join(mm.DATA_PATH, "ephemeris_files",
                           "Spitzer_ephemeris_01.dat")
data_1 = mm.MulensData(file_name=file_name_1)
data_2 = mm.MulensData(file_name=file_name_2, ephemerides_file=ephemeris_2)

# Prepare Model and Event instances
my_model = mm.Model({
    't_0': 2456836.11, 'u_0': 0.922, 't_E': 22.87,
    'pi_E_N': -0.248, 'pi_E_E': 0.234})
my_event = mm.Event([data_1, data_2], my_model, coords="17:50:00 -29:00:05")

# All get_scaled_lc() and plot_lc() calls return data scaled to the first dataset,
# unless 'dataset' keyword is used.

# Print model light curve:
print("Magnitudes:")
print(my_event.get_scaled_lc())
print("Fluxes from 2456820 to 2456850:")
print(my_event.get_scaled_lc(t_range=[2456820., 2456850.], n_epochs=20,
      phot_fmt="flux"))

# The two code lines below require discussion:
# what is default value of phot_fmt?

# RP: I think it should "scaled" to be consistent with my_event.get_scaled_lc(),
# but on the other hand phot_fmt="flux" could mean fluxes in the scale of
# that dataset, or the reference dataset and I see logic behind both of them.
# Maybe we should add "scaled_flux" and "scaled_mag" (also include both of them
# in FitData.get_residuals() to be consistent; ohh... we should also deprecate
# "scaled" in that case).

# Event.get_scaled_lc(dataset=1, phot_fmt="scaled") --> returns model magnitudes for
#                                                dataset 1 scaled to ref_data.
# Event.get_scaled_lc(dataset=data1, phot_fmt="mag") --> returns model magnitudes for
#                                                 dataset 1 scaled to that
#                                                 dataset's fs, fb.
# The above also take into account bandpass and satellite_skycoord.

# Make plot:
my_event.plot_lc()  # This one of course has many options.
plt.show()
