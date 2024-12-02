"""
Show how to plot triple lens caustics.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


raise NotImplementedError('Triple Lenses are Not Supported')

parameters = {
    't_0': 2456789.12, 'u_0': 0.01, 't_E': 34.56, 'rho': 0.0012, 'alpha': 45.,
    's_21': 1.1, 's_31': 0.5, 'q_21': 0.2, 'q_31': 0.01234, 'psi': 90.}

model = mm.Model(parameters)
model_parameters = mm.ModelParameters(parameters)

# Basic plotting:
model.plot_trajectory(caustics=True, lenses=True)  # 'lenses" is a new keyword.
# It plots dots at positions of lenses.
plt.show()

# Plot using mm.Caustics:
keys = {'s_21', 's_31', 'psi', 'q_21', 'q_31'}  # We need only these keys
caustic_parameters = {key: parameters[key] for key in keys}
# New initialization of Caustics:
caustics_1 = mm.CausticsBinary(**caustic_parameters)
caustics_1.plot()
plt.show()

# Find positions of the components and use them to plot the caustic:
positions = model_parameters.lens_positions()  # This is new function. One can
# add parameter: epoch=2456789.01 so that orbital motion is taken into account.
print(positions)  # Prints a (3,2) numpy array
caustics_2 = mm.CausticsBinary(
    q_21=parameters['q_21'], q_31=parameters['q_31'],
    lens_positions=positions)  # This keyword should also work for binary lens.
caustics_2.plot()
plt.show()

# Change position of one of the components and then plot the caustic:
# (this is why we want two ways to init Caustics for 3L events)
positions[1, 0] += 0.1
positions[1, 1] += 0.5
caustics_3 = mm.CausticsBinary(
    q_21=parameters['q_21'], q_31=parameters['q_31'],
    lens_positions=positions)
caustics_3.plot()
plt.show()
