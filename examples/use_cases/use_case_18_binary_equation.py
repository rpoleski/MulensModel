"""
Use Case 18: Binary Magnification Equation

Change how the binary magnification is calculated based on
- the time stamp of the model
- the magnification of the model (???)
"""
import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm


# Initialize the model
t_0 = 2455747.049357
t_E = 21.6796
model = mm.Model(
    {'t_0': t_0, 'u_0': 0.00352, 't_E': t_E, 'rho': 0.001647,
     'alpha': 41.35, 's': 0.5486, 'q': 0.00532})

# times to calculate the magnification
times = np.arange(t_0 - 1., t_0 + 1., 0.001)

# Set method that is used when no other method is specified
# (default value is 'point_source'):
model.set_default_magnification_method('point_source')

# Calculate the magnification using different magnification calculations
default_magnification = model.get_magnification(times)

# Specify list that give time ranges and methods:
methods = [
    2455746., 'Quadrupole', 2455746.6, 'Hexadecapole', 2455746.7, 'VBBL',
    2455747., 'Hexadecapole', 2455747.15, 'Quadrupole', 2455748.]
model.set_magnification_methods(methods)

# And specify additional parameters needed by some of the methods:
vbbl_parameters = {
    'accuracy': 0.0005}  # This is twice better than default of 0.001.
methods_parameters = {'VBBL': vbbl_parameters}
model.set_magnification_methods_parameters(methods_parameters)

# NOT IMPLEMENTED: Set times of caustic crossings. Use different
# magnification calculation based on number of source radii from the
# crossing.

accurate_magnification = model.get_magnification(times)

# Plot the differences
plt.figure()
plt.title('Magnification Curves')
label_1 = 'Default Magnification'
label_2 = 'User Specified Equations'
plt.plot(times, default_magnification, color='black', label=label_1)
plt.plot(times, accurate_magnification, color='green', label=label_2)
plt.xlabel('Times')
plt.ylabel('Magnification')
plt.legend(loc='best')

plt.show()
plt.close()

plt.figure()
plt.title('Difference in Magnification Curves')
plt.plot(times, default_magnification - accurate_magnification, color='green')
plt.xlabel('Times')
plt.ylabel('Magnification Difference')

plt.show()
