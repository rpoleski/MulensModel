"""
Use case for binary source modeling. Three types:
1. 2 sources w/rectilinear motion
2. 1 source w/xallarap
3. 2 sources w/xallarap

#1 uses Jung et al. 2017 # AJ 153, 129 (OB160733) for the fiducial values
"""
import matplotlib.pyplot as plt

import MulensModel as mm


raise NotImplementedError('Binary sources are not supported')

# Create a binary source model.
two_source_model = mm.Model({
    't_0_1': 2457501.374, 'u_0_1': 0.015, 't_0_2': 2457507.804,
    'u_0_2': 0.365, 't_E': 14.428, 'source_flux_ratio': {'I': 1.853}})

plt.figure()
two_source_model.plot_magnification()
plt.title('2-source Magnification Curve')

plt.figure()
two_source_model.plot_caustics()  # for a point lens,
# this should be a single point
two_source_model.plot_trajectory()
plt.legend()  # should automatically label the trajectories for the 2 sources.
plt.title('2-source caustics and trajectory')

# Xallarap Model
xallarap_model = mm.Model(
    {'t_0': 2457800, 'u_0': 0.001, 't_E': 32., 'xsi': [0.1, 0.2],
     'period': 3.})
xallarap_model_2 = mm.Model(
    {'t_0': 2457800, 'u_0': 0.001, 't_E': 32., 'xsi_N': 0.1, 'xsi_E': 0.2,
     'period': 3.})

plt.figure()
xallarap_model.plot_magnification()
plt.title('Xallarap Magnification Curve')

plt.figure()
xallarap_model.plot_caustics()  # for a point lens,
# this should be a single point
xallarap_model.plot_trajectory()
plt.legend()  # should automatically label the trajectories for the 2 sources.
plt.title('Xallarap trajectory')

# Two-sources with linear orbital motion see p.5 of Jung et al. 2017
# AJ 153, 129 for parameterization.
two_source_motion = mm.Model({
    't_0': 2457800, 'u_0': 0.3, 't_E': 32, 's_source': 0.2, 'q_source': 0.5,
    'alpha_source': 0.1, 'ds_dt_source': 1.0, 'dalpha_dt_source': 0.7,
    'source_flux_ratio': {'I': 0.1, 'V': 0.15}})

plt.figure()
two_source_model.plot_magnification()
plt.title('2-source Magnification Curve w/Orb Motion')

plt.figure()
two_source_model.plot_caustics()  # for a point lens,
# this should be a single point
two_source_model.plot_trajectory()
plt.legend()  # should automatically label the trajectories for the 2 sources.
plt.title('2-source caustics and trajectory w/Orb Motion')

plt.show()
