"""
Show magnification curves in different bands for a binary source event with
finite source effects and color dependent effects.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


raise NotImplementedError(
    'limb-darkening for multiple sources not implemented.')


# Create binary source models
def generate_model(finite_source=True):
    if finite_source:
        model = mm.Model({'t_0_1': 5000., 'u_0_1': 0.05,
                          't_0_2': 5005., 'u_0_2': 0.003, 'rho_2': 0.005,
                          't_E': 30.})
        model.set_limb_coeff_gamma('I', 0.5)
        model.set_limb_coeff_gamma('V', 0.7)
        model.set_magnification_methods(
            [5004., 'finite_source_LD_Yoo04', 5006.])
    else:
        model = mm.Model({'t_0_1': 5000., 'u_0_1': 0.05,
                          't_0_2': 5005., 'u_0_2': 0.003,
                          't_E': 30.})

    # model.set_source_flux_ratio_for_band('I', 0.01)
    # model.set_source_flux_ratio_for_band('V', 0.005)

    return model


# Output the magnification at time (ref_time)
def get_mag_string(model):
    ref_time = 5005.01
    # JCY - I am not actually sure how eff_mag should work for finite
    # source models if a band is not specified. Maybe it should fail...
    # eff_mag = model.get_magnification(ref_time)
    eff_mag_i = model.get_magnification(
        ref_time, bandpass='I', source_flux_ratio=0.01)
    eff_mag_v = model.get_magnification(
        ref_time, bandpass='V', source_flux_ratio=0.005)
    (mag_1_i, mag_2_i) = model.get_magnification(
        ref_time, bandpass='I', separate=True)
    (mag_1_v, mag_2_v) = model.get_magnification(
        ref_time, bandpass='V', separate=True)
    str_1 = '{0:8.4f} {1:8.4f} '.format(
        eff_mag_i, eff_mag_v)
    str_2 = '{0:8.4f} {1:8.4f} {2:8.4f} {3:8.4f}'.format(
        mag_1_i, mag_2_i, mag_1_v, mag_2_v)
    string = str_1 + str_2

    return string


model_fs = generate_model(finite_source=True)
model_ps = generate_model(finite_source=False)

print('{0:6} {1:8} {2:8} {3:8} {4:8} {5:8} {6:8} {7:8}'.format(
    'model', 'eff_mag', 'eff_mag_i', 'eff_mag_v', 'mag_1_i', 'mag_2_i',
    'mag_1_v', 'mag_2_v'))
print('{0:6} {1}'.format('FS', get_mag_string(model_fs)))
print('{0:6} {1}'.format('PS', get_mag_string(model_ps)))

# Show the different light curves
plt.figure()
plt.title('With finite source effects')
model_fs.plot_magnification(
    bandpass='I', color='red', linestyle='-', label='I, FS')
model_fs.plot_magnification(
    bandpass='V', color='blue',  linestyle='-', label='V, FS')
model_ps.plot_magnification(
    bandpass='I', color='red', linestyle=':', label='I, PS')
model_ps.plot_magnification(
    bandpass='V', color='blue',  linestyle=':', label='V, PS')
plt.legend(loc='best')
plt.show()
