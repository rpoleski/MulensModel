"""
Example use of MulensModel to model a single source binary lens with an
external mass sheet, as in Peirson et al. 2022 ApJ 927, 24
and Vedantham et al. 2017 ApJ 845, 89.
Binary lens with external mass sheet assumes a point source when using VBBL.
"""

import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm


# Define lens model and source parameters
t_0 = 300
u_0 = 0.25
t_E = 500
s = 1.0
q = 0.4
alpha = 0
K = 0.0
G = complex(0.0, -0.2)

time = np.arange(t_0-775., t_0+775., 0.5, dtype=float)

lens = mm.model.Model({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,
        'convergence_K': K, 'shear_G': G})
no_shear = mm.model.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q,
                           'alpha': alpha})

# Plot magnification curve and caustics
(fig, (ax1, ax2)) = plt.subplots(figsize=(10, 5), ncols=2)

ax1.plot(time, lens.get_magnification(time), color='r')
ax1.plot(time, no_shear.get_magnification(time), alpha=0.4)

ax2.set_xlim(-1.5, 1.5)
ax2.set_ylim(-1.5, 1.5)
ax2.set_xlabel("x", fontweight="bold")
ax2.set_ylabel("y", fontweight="bold")
no_shear.plot_trajectory(t_range=[t_0 - 475, t_0], caustics=True,
                         color='green')
mm.Caustics(s=s, q=q).plot(alpha=0.3)
lens.plot_trajectory(t_range=[t_0 - 475, t_0], caustics=True, color='green')
lens.plot_trajectory(t_range=[t_0, t_0 + 475], color='blue')

plt.show()
