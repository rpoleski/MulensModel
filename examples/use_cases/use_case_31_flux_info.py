"""
use_case_31_flux_info.py

Adapted from use_case_25_flux_constraint.py

Shows different ways to access fitted fluxes.
"""
import os
import MulensModel
import numpy as np


# Functioniality to set the default zeropoint to 18
# (more appropriate for ground-based data.)
MulensModel.utils.MAG_ZEROPOINT = 18.

# Add data from OB161195
datasets = []
file_names = ['KCT01I.dat', 'KCT41I.dat', 'KCT42I.dat', 'KSA01I.dat',
              'KSA41I.dat', 'KSA42I.dat', 'KSS01I.dat', 'KSS41I.dat',
              'KSS42I.dat', 'spitzer_b12.dat']
dir_ = os.path.join(MulensModel.DATA_PATH, "photometry_files", "OB161195")
for file_name in file_names:
    file_ = os.path.join(dir_, file_name)
    datasets.append(MulensModel.MulensData(
        file_name=file_, add_2450000=True,
        plot_properties={'label': file_name}))

# Close-- model
model = MulensModel.Model(
    {'t_0': 2457568.7692, 'u_0': -0.05321, 't_E': 9.96, 'rho': 0.00290,
     'pi_E_N': -0.2154, 'pi_E_E': -0.380,
     'alpha': np.rad2deg(-0.9684), 's': 0.9842, 'q': 0.0000543})

methods = [2457560., 'VBBL', 2457580.]

model.set_magnification_methods(methods)
model.set_default_magnification_method('point_source_point_lens')

event = MulensModel.Event(datasets=datasets, model=model,
                          coords="17:55:23.50 -30:12:26.1")

print(event.get_ref_fluxes())

# Output the source and blend fluxes in a "nice" way.
print('Observatory, Source Flux, Blend Flux:')
for dataset in event.datasets:
    (source_flux, blend_flux) = event.get_flux_for_dataset(dataset)
    print('{0:15} {1:8.2f} {2:8.2f}'.format(
        dataset.plot_properties['label'], source_flux[0], blend_flux))
