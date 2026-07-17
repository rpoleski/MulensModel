"""
Here we show how to define basic xallarap models.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


def make_plots(model, title, source_flux_ratio=None):
    """
    Make plots of magnification and trajectory.
    For models with two luminous sources, one has to provide source_flux_ratio. 
    """
    model.plot_magnification(source_flux_ratio=source_flux_ratio)
    plt.title(title + " - magnification")

    plt.figure()
    model.plot_trajectory(caustics=True)
    plt.title(title + " - trajectory")
    plt.axis('equal')

    plt.show()


# Define set of parameters.
# First three basic ones: t_0, u_0, t_E.
# Then five parameters that define the circular xallarap orbit.
# I'm chosing xi_semimajor_axis comparable to u_0 and xi_period shorter than t_E,
# so that the xallarap effect dominates microlensing signal.
parameters = {
    't_0': 5002., 'u_0': 0.4567, 't_E': 35.,
    'xi_period': 12.345, 'xi_semimajor_axis': 0.25, 'xi_Omega_node': 0.123,
    'xi_inclination': 0., 'xi_argument_of_latitude_reference': 24.68}

# Make a model and plot it:
model_1 = mm.Model(parameters)
make_plots(model_1, 'circular 1S1L')

# Change that into an eccentric orbit - it needs two additional parameters (eccentricity and orbit orientation).
parameters_ecc = {**parameters, 'xi_eccentricity': 0.8, 'xi_omega_periapsis': 180.}
model_2 = mm.Model(parameters_ecc)
make_plots(model_2, 'eccentric 1S1L')

# Once more make model with circular orbit, but this time both components are bright.
# First we need to add source mass ratio.
model_3 = mm.Model({**parameters, 'q_source': 0.3})
# For plotting we also have to provide flux ratio of source, so that the effective magnification can be calculated.
make_plots(model_3, 'circular 2S1L', 0.2)

# Final model - eccentric orbit and two sources.
model_4 = mm.Model({**parameters_ecc, 'q_source': 0.3})
make_plots(model_4, 'eccentric 2S1L', 0.2)

