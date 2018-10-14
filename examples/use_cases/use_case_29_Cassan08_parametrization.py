"""
The paper:

Cassan 2008 http://adsabs.harvard.edu/abs/2008A%26A...491..587C
"An alternative parameterisation for binary-lens caustic-crossing events"

introduced a different way to parametrize caustic-crossing binary lens models.
Instead of t_0, u_0, t_E, and alpha parameters, the new parameterization
uses two epochs of caustic crossing and corresponding two curvilinear
abscissa along the caustic. The epochs of caustic crossing are well defined
in typical cases and it is advantageous to use them as fitting parameters.

These curvilinear abscissa are normalized so that XXX
"""
import matplotlib.pyplot as plt

import MulensModel as MM


model = MM.Model({'s': 1.01, 'q': 0.001, 'rho': 0.002,
                  's_caustic_in': 1.36, 's_caustic_out': 0.72,
                  't_caustic_in': 6543.123, 't_caustic_out': 6550.987}
# For models with orbital motion we will be using also:
# s_in s_out
# which are simply separation values at t_caustic_in and t_caustic_out.
# The values of s_in and s_out will be not much different from s.

print(model.n_lenses, model.n_sources) # Returns "2 1"

# Plot trajectory and caustic:
model.plot_trajectory(caustics=True)
plt.show()

print("The same model in standard parametrization:")
print("t_0 :", model.parameters.t_0)
print("u_0 :", model.parameters.u_0)
print("t_E :", model.parameters.t_E)
print("alpha :", model.parameters.alpha)
print("Other parameters are the same.")

print("You may also be interested in t_star:")
print(model.parameters.t_star) 
# This makes sure that internally t_E is properly defined.

