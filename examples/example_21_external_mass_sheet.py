"""
Example use of MulensModel to model a single source binary lens with an
external mass sheet, as in Peirson et al. 2022 ApJ 927, 24
and Vedantham et al. 2017 ApJ 845, 89.
Binary lens with external mass sheet assumes a point source when using VBBL.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


standard_parameters = {
    't_0': 300., 'u_0': 0.25, 't_E': 500, 's': 1., 'q': 0.4, 'alpha': 180.}
mass_sheet = {'convergence_K': 0.0, 'shear_G': complex(0.0, -0.2)}

shear_model = mm.model.Model({**mass_sheet, **standard_parameters})
no_shear_model = mm.model.Model(standard_parameters)

(_, (ax1, ax2)) = plt.subplots(figsize=(10, 5), ncols=2)
t_range = [standard_parameters['t_0'] - 300, standard_parameters['t_0'] + 550]

# Magnification plot:
plt.sca(ax1)
shear_model.plot_magnification(t_range=t_range, color='r')
no_shear_model.plot_magnification(t_range=t_range, color='b')

# Caustics and trajectory plot:
plt.sca(ax2)
shear_model.plot_trajectory(t_range=t_range, caustics=True, color='green')
no_shear_model.plot_caustics(color='b')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

plt.show()
