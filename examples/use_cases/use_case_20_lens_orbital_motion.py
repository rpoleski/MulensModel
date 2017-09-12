from astropy import units as u
import matplotlib.pyplot as plt

import MulensModel

"""
Use case presenting binary lens orbital motion models. Also comparison with static 
binary model is provided.
"""

# point lens parameters:
t_0 = 2457123.456
u_0 = 0.0345
t_E = 30.00 * u.day
rho = 0.00345

# binary lens parameters:
q = 0.123
alpha_0 = 12.345 * u.deg
dalpha_dt = 50. u.deg / u.year
s_0 = 1.5
ds_dt = 0.5 / u.year

#Generate a model.
model = MulensModel.Model()
model.set_parameters(t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, q=q, 
                    alpha_0=alpha_0, dalpha_dt=dalpha_dt, 
                    s_0=s_0, ds_dt=ds_dt) 
                    # t_0_kep is not provided hence defaults to t_0

model_static = MulensModel.Model()
model_static.set_parameters(t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, q=q,
                    alpha=alpha_0, s=s_0)

dt = 36.525 # This is in days.

# Get the values of parameters:
# print(model.s) - this would raise an exception.
print(model.s_0) 
print(model.s_for_epoch(t_0)) # Prints the same as previous one.
print(model.s_for_epoch(t_0+dt)) # should return 1.55
print(model.alpha_for_epoch(epoch=t_0-dt)) # should return 7.345 u.deg
# for static model all three commands give the same value:
print(model_static.s)
print(model_static.s_for_epoch(t_0)) # Yes, s_for_epoch() works for both static and orbiting models.
print(model_static.s_for_epoch(t_0+dt))
# print(model_static.s_0) - this would raise an exception.

# Print projected orbital velocity
print(model.gamma_parallel) # should return 0.3333333 1/u.year
print(model.gamma_perp) # should return -50 1/u.year or -0.87266 u.rad/u.year 
# (the minus sign comes from the definition in Skowron et al. 2011)
print(model.gamma) # should return 0.9346 1/u.year

# Make a nice plot
plt.figure()
model.plot_caustics(epoch=t_0)
model.plot_caustics(epoch=t_0+dt, c='g') # second caustics are green
model.plot_trajectory()
t = 'This plot shows a nice curved trajectory and caustics for 2 different epochs'
plt.title(t)
plt.show()

