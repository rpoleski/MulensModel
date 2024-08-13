"""
Show several examples of many-body models

Parameters correspond to OB151459: Hwang et al. 2018, AJ, 155, 259
"""
import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np

# 1L3S event
subhead_1L3S = 'Table 3: Static 1L3S Model '

params_1L3S = {
    't_E': 4.921,
    't_0_1': 7199.946, 'u_0_1': 0.065,
    't_0_2': 7200.193, 'u_0_2': 2.638e-3, 'rho_2': 4.503e-3,
    't_0_3': 7200.202, 'u_0_3': 0.281e-3, 'rho_3': 0.631e-3}

flux_ratios_1L3S = {'q_2': 0.014, 'q_3': 0.006}
source_flux_1_1L3S = 0.056
source_flux_2_1L3S = source_flux_1_1L3S * flux_ratios_1L3S['q_2']
source_flux_3_1L3S = source_flux_1_1L3S * flux_ratios_1L3S['q_3']
blend_flux_1L3S = 0.100
source_fluxes_1L3S = [source_flux_1_1L3S, source_flux_2_1L3S, source_flux_3_1L3S]

mag_methods_1L3S = [7200.1, 'finite_source_uniform_Gould94', 7200.3]
model_1L3S = mm.Model(params_1L3S)
model_1L3S.set_magnification_methods(mag_methods_1L3S, 2)
model_1L3S.set_magnification_methods(mag_methods_1L3S, 3)

# 2L2S event
subhead_2L2S = 'Table 5: 2L2S Model - close '

t_E = 4.613
s = 0.857
q = 2.374e-3
alpha = -np.rad2deg(0.890)
params_2L2S = {
    't_E': t_E, 's': s, 'q': q, 'alpha': alpha,
    't_0_1': 7199.941, 'u_0_1': 0.069, 'rho_1': 0.688e-3,
    't_0_2': 7200.057, 'u_0_2': 0.039, 'rho_2': 0.542e-3}

q_F_2L2S = 0.181
source_flux_1_2L2S = 0.061
source_flux_2_2L2S = source_flux_1_2L2S * q_F_2L2S
blend_flux_2L2S = 0.095
source_fluxes_2L2S = [source_flux_1_2L2S, source_flux_2_2L2S]

mag_methods_2L2S = [7200.1, 'VBBL', 7200.3]

model_2L2S = mm.Model(params_2L2S)
model_2L2S.set_magnification_methods(mag_methods_2L2S)

# Make the plots
times_perturb_1 = np.arange(7200.1, 7200.3, 0.001)
times_perturb_2 = np.arange(7200.195, 7200.21, 0.0001)
times_lc = np.arange(7190., 7210., 0.01)
times = np.hstack((times_lc, times_perturb_1, times_perturb_2))
times.sort()

plt.figure(figsize=(8, 4))
plt.suptitle('OB151459')

plt.subplot(1, 2, 1)
plt.title(subhead_1L3S)
model_1L3S.plot_trajectory()
plt.gca().set_aspect('equal')
plt.xlim(-0.1, 0.1)
plt.ylim(-0.1, 0.1)
plt.minorticks_on()

plt.subplot(1, 2, 2)
plt.title(subhead_2L2S)
model_2L2S.plot_trajectory(caustics=True)
plt.gca().set_aspect('equal')
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.minorticks_on()

plt.figure()
plt.title('OB151459')
model_1L3S.plot_lc(times, source_flux=source_fluxes_1L3S, blend_flux=blend_flux_1L3S, label=subhead_1L3S, color='black')
model_2L2S.plot_lc(times, source_flux=source_fluxes_2L2S, blend_flux=blend_flux_2L2S, label=subhead_2L2S, color='red')
plt.legend()
plt.minorticks_on()
plt.tight_layout()

plt.show()
