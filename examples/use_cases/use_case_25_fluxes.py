"""
Use cases showing how to implement a flux constraint.
"""
import numpy as np
import glob

import MulensModel


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

    f_s_ogle = event.get_source_flux(dataset=KMT)
    f_s_spitz = event.get_source_flux(dataset=Spitzer) 

    color = -2.5 * np.log10(f_s_ogle / f_s_spitz)
    
    return (color - source_color)**2 / sigma**2

# Add data from OB161195
datasets = []
file_names = ['KCT01I.dat', 'KCT41I.dat', 'KCT42I.dat', 'KSA01I.dat',
            'KSA41I.dat', 'KSA42I.dat', 'KSS01I.dat', 'KSS41I.dat', 
            'KSS42I.dat', 'spitzer_b12.dat']
for file_ in file_names:
    datasets.append(MulensModel.MulensData(
            file_name=os.path.join(
        MulensModel.MODULE_PATH, "data", "photometry_files", "OB161195", file_)))

# Close-- model
model = MulensModel.Model(
    {'t_0': 2457568.7692, 'u_0': -0.05321, 't_E': 9.96, 'rho': 0.00290, 
     'pi_E_N': -0.2154, 'pi_E_E': -0.380, 
     'alpha': np.rad2deg(-0.9684), 's': 0.9842, 'q': 0.0000543})

event = MulensModel.Event(datasets=datasets, model=model)
print(chi2_fun(event))

