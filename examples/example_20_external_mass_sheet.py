"""
Example use of MulensModel to model a single source binary lens with an
external mass sheet, as in Peirson et al. 2022 ApJ, 
Vedantham et al. 2017 ApJ 845, 89. 
Binary lens with external mass sheet assumes a point source when using VBBL.
"""

import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm

#Define lens model and source parameters
s = 1.0
q = 0.01
alpha = 270 
K = 0.1
G = complex(0.05,0.0)
t_0 = 300
t_E = 500
rho = 1e-5
u_0 = -0.1

Day = np.arange(-500,1000,1,dtype=float)

lens = mm.model.Model(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha, 'convergence_K': K, 
        'shear_G': G, 'rho':rho})

lens.set_magnification_methods([min(Day), 'vbbl', max(Day)])

#Plot magnification curve and caustics
fig, (ax1,ax2) = plt.subplots(figsize=(10,5),ncols=2)

ax1.plot(Day, lens.get_magnification(Day),color='r')

ax2.set_xlim(-1,1)
ax2.set_ylim(-1,1)
ax2.set_xlabel("x", fontweight="bold")
ax2.set_ylabel("y", fontweight="bold")
lens.plot_trajectory(t_range=[t_0 - t_E, t_0], caustics=True, color='red',)
lens.plot_trajectory(t_range=[t_0, t_0 + t_E], caustics=True, color='blue',)