"""
use_case_31_flux_info.py

Adapted from use_case_25_flux_constraint.py
New functionality marked by # *** NEW ***

Shows different ways to access fitted fluxes.
"""
import os
import MulensModel

raise NotImplementedError('This use case has not been implemented.')

# Function to set the default zeropoint to 18
# (more appropriate for ground-based data.)
MulensModel.set_mag_zeropoint(18.)

# Add data from OB161195
datasets = []
file_names = ['KCT01I.dat', 'KCT41I.dat', 'KCT42I.dat', 'KSA01I.dat',
              'KSA41I.dat', 'KSA42I.dat', 'KSS01I.dat', 'KSS41I.dat',
              'KSS42I.dat', 'spitzer_b12.dat']
dir_ = os.path.join(MulensModel.DATA_PATH, "photometry_files", "OB161195")
for file_name in file_names:
    file_ = os.path.join(dir_, file_name)
    datasets.append(MulensModel.MulensData(
        file_name=file_, add_2450000=True, plot_properties={'label': file_name}))

# Close-- model
model = MulensModel.Model(
    {'t_0': 2457568.7692, 'u_0': -0.05321, 't_E': 9.96, 'rho': 0.00290,
     'pi_E_N': -0.2154, 'pi_E_E': -0.380,
     'alpha': np.rad2deg(-0.9684), 's': 0.9842, 'q': 0.0000543})

methods = [7560., 'VBBL', 7580.]

model.set_magnification_methods(methods)
model.set_default_magnification_method('point_source_point_lens')

event = MulensModel.Event(datasets=datasets, model=model,
                          coords="17:55:23.50 -30:12:26.1")

# *** NEW ***
print(event.get_ref_fluxes())
print(event.fluxes)
# should be a (n, m) array where n = the number of datasets, and m = n_sources + 1
# --> one flux value for each source + the blend flux. It makes sense for the blend
# flux to be last because then it can be consistently accessed with index = -1.
print(event.source_fluxes)
print(event.blend_fluxes)

# Output the source and blend fluxes in a "nice" way.
print('Observatory, Source Flux, Blend Flux')
for i in range(len(datasets)):
    print('{0:10} {1:8.2f} {2:8.2f}'.format(
        event.datasets[i].plot_properties['label'], event.fits[ datasets[i] ].source_flux,
        event.fits[ datasets[i] ].blend_flux))

# *** END NEW ***