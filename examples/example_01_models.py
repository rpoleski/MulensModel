from MulensModel.model import Model

import matplotlib.pyplot as pl
import numpy as np

# Create a PSPL model
t_0 = 3583.
u_0 = 0.3
t_E = 12.

pspl = Model(t_0=t_0, u_0=u_0, t_E=t_E)

#Create a planet model with same PSPL parameters
s = 1.5
q = 0.001
alpha = np.rad2deg(0.37) 
print('WARNING: This value of alpha does not produce a planet perturbation. JCY.')

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

#Plot planet and PSPL models and show difference in magnification at planet perturbation
pl.figure()
pspl.plot_magnification(color='blue', linestyle=':', zorder=1)
planet.plot_magnification(color='red',linestyle='-', zorder=2)
t_p = t_0+0.768*t_E
pl.plot(
    [t_p, t_p], [pspl.magnification(t_p), planet.magnification(t_p)], 
    color='black')
pl.title('Planet vs. Point Lens Models')

pl.show()
