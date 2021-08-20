"""
Plot source positions.
"""
from matplotlib import pyplot as plt
import numpy as np

import MulensModel as mm


# Define Model object
params = {"t_0": 0.0, "u_0": 0.1, "t_E": 25.0, "rho": 1e-2, "s": 1.0,
          "q": 1e-3, "alpha": 90.}
model = mm.Model(parameters=params)

# Create figure environment
fig = plt.figure(figsize=(10, 8), constrained_layout=False)
ax = fig.add_subplot(111)

# Plot caustics and source trajectory for specified model
model.plot_caustics(color="black")
trajectory_kwargs = {'t_start': params["t_0"]-0.15*params["t_E"],
                     't_stop': params["t_0"]+0.15*params["t_E"],
                     'caustics': False, 'color': "blue"}
model.plot_trajectory(**trajectory_kwargs)

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
plt.axis('equal')  # So that circles don't look like ellipses.
plt.show()

# Repeat above, but also plot data:
epochs = np.array([-3.5, -1., 0.5, 0.6, 0.75, 0.95, 2.])
dataset = mm.MulensData([epochs, epochs*0.+20., epochs*0.+.1])
event = mm.Event(datasets=dataset, model=model)
event.plot_trajectory(**trajectory_kwargs)
event.model.plot_caustics(color='black')
event.model.plot_source(times, **kwargs)
event.plot_source_for_datasets()
plt.axis('equal')
plt.show()

# Plotting for binary source model:
model = mm.Model({
    't_0_1': 5000., 'u_0_1': 0.005, 'rho_1': 0.001,
    't_0_2': 5030., 'u_0_2': 0.0003, 't_star_2': 0.03, 't_E': 25.})
times = np.linspace(4980., 5050.)
model.plot_trajectory()
model.plot_source(times)
plt.show()  # Circles are squished because of X-axis range.

# Same as above, but no source size provided:
model = mm.Model({
    't_0_1': 5000., 'u_0_1': 0.005,
    't_0_2': 5030., 'u_0_2': 0.0003, 't_E': 25.})
model.plot_trajectory()
model.plot_source(times)
plt.show()
