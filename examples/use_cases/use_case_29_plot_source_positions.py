from matplotlib import pyplot as plt
import numpy as np

import MulensModel as MM


# Define Model object
params = {"t_0": 0.0, "u_0": 0.1, "t_E": 25.0, "rho": 1e-2, "s": 1.0,
          "q": 1e-3, "alpha": 90.}
model = MM.Model(parameters=params)

# Create figure environment
fig = plt.figure(figsize=(10, 8), constrained_layout=False)
ax = fig.add_subplot(111)

# Plot caustics and source trajectory for specified model
model.plot_caustics(color="black")
model.plot_trajectory(
    t_start=params["t_0"]-0.15*params["t_E"],
    t_stop=params["t_0"]+0.15*params["t_E"],
    caustics=False, color="blue")

# Stages for plotting source positions (MM documentation has this listed but
# describes it as not yet implemented)
# (1) plot source position along trajectory for specified time(s)
# (2) plot source position at arbitrary (x, y) [theta_E] position(s)
# (3) plot source positions along trajectory, color-coded to time(s) of
#     dataset(s) (requires MulensData)
# ...allow for color and marker style (open circle; maybe "x" as well?)
times = params["t_0"] - 0.05 * params["t_E"]
kwargs = {}  # You can add some kwargs here and they will be passed
# to plt function.
model.plot_source(times, **kwargs)
plt.axis('equal') # So that circles don't look like elipses.
plt.show()

# Plotting for binary source model:
model = MM.Model({
    't_0_1': 5000., 'u_0_1': 0.005, 'rho_1': 0.001,
    't_0_2': 5030., 'u_0_2': 0.0003, 't_star_2': 0.03, 't_E': 25.})
times = np.linspace(4980., 5050.)
model.plot_trajectory()
model.plot_source(times)
plt.show()

# Same as above, but no source size provided:
model = MM.Model({
    't_0_1': 5000., 'u_0_1': 0.005,
    't_0_2': 5030., 'u_0_2': 0.0003, 't_E': 25.})
model.plot_trajectory()
model.plot_source(times)
plt.show()

