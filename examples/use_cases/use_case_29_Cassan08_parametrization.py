"""
The paper:

Cassan 2008 http://adsabs.harvard.edu/abs/2008A%26A...491..587C
"An alternative parameterisation for binary-lens caustic-crossing events"

introduced a different way to parametrize caustic-crossing binary lens models.
Instead of t_0, u_0, t_E, and alpha parameters, the new parameterization
uses two epochs of caustic crossing and corresponding two curvilinear
abscissa along the caustic. The epochs of caustic crossing are well defined
in typical cases and it is advantageous to use them as fitting parameters.
"""
import matplotlib.pyplot as plt

import MulensModel as MM


model = MM.Model({'s': 1.01, 'q': 0.001, 'rho': 0.002,
                  'x_caustic_in': 0.36, 'x_caustic_out': 0.72,
                  't_caustic_in': 6543.123, 't_caustic_out': 6550.987})
# For models with orbital motion we will be using also:
# s_in s_out
# which are simply separation values at t_caustic_in and t_caustic_out.
# The values of s_in and s_out will be not much different from s.

print(model.n_lenses, model.n_sources) # Returns "2 1".

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

# The function below only runs calculations that speed-up further calculations.
# It significantly speeds-up model fitting but has no immediate result.
# Maybe it's possible to run these calculations once and save everything in
# a file? Then this function won't be needed.
model.precompute_Cassan08_normalization(
    x_caustic_start=0.9, x_caustic_stop=1.1, q_start=1.e-5, q_stop=0.1)
# There may be additional parameters that specify how finely the values are
# calculated.

# The function below calculates the integral for selected (s, q) and
# remembers it, so that further calculations of (t_0, u_0, t_E, alpha)
# with these (s, q) would require just interpolation instead of integration.
model.precompute_Cassan08_for_fixed_s_q(s=1.01, q=0.001)

