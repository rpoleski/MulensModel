import matplotlib.pyplot as pl
import numpy as np

import MulensModel
from MulensModel.model import Model

# Create a PSPL model
t_0 = 3583.
u_0 = 0.3
t_E = 12.

pspl = Model(t_0=t_0, u_0=u_0, t_E=t_E)

#Create a planet model with same PSPL parameters
s = 1.5
q = 0.001
alpha = np.rad2deg(np.pi - 0.37)

planet = Model(t_0=t_0, u_0=u_0, t_E=t_E, s=s, q=q, alpha=alpha)

#Plot PSPL model
pl.figure()
pspl.plot_magnification()
pl.title('A Point Lens Model')

#Plot PSPL model in magnitudes with arbitrary blending
pl.figure()
pspl.plot_lc(f_source=1.0, f_blend=0.0, label='fs=1.0, fb=0.0')
pspl.plot_lc(f_source=0.5, f_blend=0.5, label='fs=0.5, fb=0.5')
pl.legend(loc='best')
pl.title('A Point Lens Model in Magnitudes')

#Plot planet and PSPL models and show difference in magnification at
#planet perturbation
pl.figure()
pspl.plot_magnification(
    color='blue', linestyle=':', zorder=1, label='Point Lens')
planet.plot_magnification(
    color='red',linestyle='-', zorder=2, label='Planet')
t_p = t_0+0.768*t_E
pl.title('Planet vs. Point Lens Models')
pl.legend(loc='best')

#Plot detail of the planet perturbation
pl.figure()
planet.plot_magnification(
    t_range=[3592, 3593], 
    color='red',linestyle='-', zorder=2, label='Planet')
pl.title('Planetary Perturbation Detail')

pl.show()
