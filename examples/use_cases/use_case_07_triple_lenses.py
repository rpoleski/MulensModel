"""
Triple lens definition.
"""
from matplotlib import pyplot as plt

import MulensModel as mm


raise NotImplementedError('Triple Lenses are Not Supported')

model = mm.Model({
    't_0': 2456789.12, 'u_0': 0.01, 't_E': 34.56, 'rho': 0.0012, 'alpha': 45.,
    's_21': 1.1, 's_31': 0.5, 'q_21': 0.2, 'q_31': 0.01234, 'psi': 90.})
# The line above is shows the names of triple lens parameters.

model.set_default_magnification_method('point_source_point_lens')
model.set_magnification_methods([
    2456700., 'point_source', 2456730., 'hexadecapole', 2456800.])

model.plot_magnification(subtract_2450000=True)
plt.show()

model.plot_trajectory(caustics=True)
plt.show()

# Previously we wanted to have Model.parameters.s which should return a list
# (e.g., [1.2, 0.5]) but it doesn't seem to be a good idea now.
