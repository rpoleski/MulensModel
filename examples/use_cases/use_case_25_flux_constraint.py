"""
Use cases showing how to implement a flux constraint.

Shows how to access fitted fluxes using event.source_fluxes
in contrast to use_case_32, which shows how to access fitted fluxes
using event.fits.
"""
import numpy as np
import os

import MulensModel as mm


def chi2_fun(theta, event, parameters_to_fit):
    """
    for given event set attributes from parameters_to_fit (list of
    str) to values from theta list
    """
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])

    chi2 = event.get_chi2()
    constraint = get_color_constraint(event)

    return chi2 + constraint


def get_color_constraint(event):
    """
    Calculate the color constraint chi2 penalty

    KMT = *int*, dataset number for KMTC I-band
    Spitzer = *int*, dataset number for Spitzer
    """
    KMT = 0
    Spitzer = 9

    # Color constraint for OB161195 (I_KMT - L_Spitzer)
    (source_color, sigma) = (0.78, 0.03)

    (f_s_kmt, f_b_kmt) = event.get_flux_for_dataset(datasets[KMT])
    (f_s_spitzer, f_b_spitzer) = event.get_flux_for_dataset(datasets[Spitzer])

    color = -2.5 * np.log10(f_s_kmt[0] / f_s_spitzer[0])

    return (color - source_color)**2 / sigma**2


# Add data from OB161195
datasets = []
file_names = ['KCT01I.dat', 'KCT41I.dat', 'KCT42I.dat', 'KSA01I.dat',
              'KSA41I.dat', 'KSA42I.dat', 'KSS01I.dat', 'KSS41I.dat',
              'KSS42I.dat', 'spitzer_b12.dat']
dir_ = os.path.join(mm.DATA_PATH, "photometry_files", "OB161195")
spitzer_ephemeris = os.path.join(
            mm.DATA_PATH, "ephemeris_files", "Spitzer_ephemeris_02.dat")
for file_name in file_names:
    file_ = os.path.join(dir_, file_name)
    if file_name[0:7] == "spitzer":
        ephemerides_file = spitzer_ephemeris
    else:
        ephemerides_file = None

    datasets.append(mm.MulensData(
        file_name=file_, add_2450000=True, ephemerides_file=ephemerides_file))

# Close-- model
model = mm.Model(
    {'t_0': 2457568.7692, 'u_0': 0.05321, 't_E': 9.96, 'rho': 0.00290,
     'pi_E_N': -0.2154, 'pi_E_E': -0.380,
     'alpha': np.rad2deg(-0.9684), 's': 0.9842, 'q': 0.0000543})

methods = [2457560., 'VBBL', 2457580.]

model.set_magnification_methods(methods)
model.set_default_magnification_method('point_source_point_lens')

coords = "17:55:23.50 -30:12:26.1"
event = mm.Event(datasets=datasets, model=model, coords=coords)

print(chi2_fun([], event, []))
