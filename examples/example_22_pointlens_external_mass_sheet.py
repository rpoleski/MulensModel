"""
Example use of MulensModel to model a PSPL with an
external mass sheet, Chang-refsdal.
PSPL with external mass sheet assumes a point source.
"""
import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm


# Define lens model and source parameters
alpha = 0
K = 0.1
G = complex(0.1, -0.05)
t_0 = 300
t_E = 500
u_0 = 0.1
d_t = 475.

time = np.arange(t_0-d_t, t_0+d_t, 0.5, dtype=float)

no_shear = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
with_shear = mm.Model({**no_shear.parameters.parameters,
                       'convergence_K': K, 'shear_G': G, 'alpha': alpha})

(_, (ax1, ax2)) = plt.subplots(figsize=(10, 5), ncols=2)

# Plot magnification curves:
plt.sca(ax1)
with_shear.plot_magnification(times=time, color='r')
no_shear.plot_magnification(times=time, alpha=0.4)

# Plot trajectories and a caustic:
plt.sca(ax2)
with_shear.plot_trajectory(t_range=[t_0 - d_t, t_0], caustics=True,
                           color='green')
with_shear.plot_trajectory(t_range=[t_0, t_0 + d_t], color='blue')
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')

plt.show()
