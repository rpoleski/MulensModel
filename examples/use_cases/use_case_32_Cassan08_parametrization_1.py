"""
The paper:

Cassan 2008 https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract
"An alternative parameterisation for binary-lens caustic-crossing events"

introduced a different way to parametrize caustic-crossing binary lens models.
Instead of t_0, u_0, t_E, and alpha parameters, the new parameterization
uses two epochs of caustic crossing and corresponding two curvilinear
abscissa along the caustic. The epochs of caustic crossing are well defined
in typical cases and it is advantageous to use them as fitting parameters.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


model = mm.Model({'s': 1.5, 'q': 0.2, 'rho': 0.00271,
                  'x_caustic_in': 0.16698, 'x_caustic_out': 0.81166,
                  't_caustic_in': 6021.43701, 't_caustic_out': 6056.25338})

print("Number of lenses and sources:")
print(model.n_lenses, model.n_sources)  # Returns "2 1".

print("\nThe same model in standard parametrization:")
print("t_0 :", model.parameters.t_0)
print("u_0 :", model.parameters.u_0)
print("t_E :", model.parameters.t_E)
print("alpha :", model.parameters.alpha)
print("Other parameters are the same.")

print("You may also be interested in t_star:")
print(model.parameters.t_star)
# This makes sure that internally t_E is properly defined.

# Plot trajectory and caustic:
model.plot_trajectory(caustics=True)
plt.show()
# and the magnification curve:
model.set_magnification_methods([5900., 'VBBL', 6200.])
model.plot_magnification(dt=0.01)
plt.show()
