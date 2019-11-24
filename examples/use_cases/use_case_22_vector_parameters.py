"""
Model parameters accessed as a vector
"""
import MulensModel as mm


raise NotImplementedError('ModelParameters.vector not implemented')

# JCY: I don't see the point of returning a vector. I also think it
# creates a risk for the user: there is no way to identify which value
# goes with which variable, you have to know this in advance.

# Show a new way of accessing model parameters.
t_0 = 2456789.0123
u_0 = 0.1
t_E = 15.
rho = 0.001
pi_E_N = 0.1
pi_E_E = 0.2

parameters = mm.ModelParameters(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E,
     'pi_E_N': pi_E_N, 'pi_E_E': pi_E_E})
print(parameters.vector)
# returns np.array([t_0, u_0, t_E, pi_E_N, pi_E_E])

# Below is non-functional:
# A new variable cannot be introduced after the model has been defined.
#
# JCY: This exactly illustrates the confusion I'm worried
# about. Why is rho introduced in the middle of the vector? This is an
# arbitrary convention commonly used by microlensers, but there is no
# fundamental reason why rho should come after t_E rather than after
# parallax. It is more transparent just to print (and access) the dictionary.
parameters.rho = rho
print(parameters.vector)
# returns np.array([t_0, u_0, t_E, rho, pi_E_N, pi_E_E])
