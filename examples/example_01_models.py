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
pspl.plot()
pl.title('A Point Lens Model')

#Plot planet and PSPL models and show difference in magnification at planet perturbation
pl.figure()
pspl.plot(color='blue', linestyle=':', zorder=1)
planet.plot(color='red',linestyle='-', zorder=2)
t_p = t_0+0.768*t_E
pl.plot([t_p,t_p], [pspl.magnification(t_p),planet.magnification(t_p)], color='black')
pl.title('Planet vs. Point Lens Models')

pl.show()
