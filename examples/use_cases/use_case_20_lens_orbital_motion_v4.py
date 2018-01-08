from astropy import units as u
import matplotlib.pyplot as plt

import MulensModel

"""
Use case presenting binary lens orbital motion models. Also comparison
with static binary model is provided.
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
model_orb = MulensModel.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho, 'q': q, 
     'alpha': alpha_0, 'dalpha_dt': dalpha_dt, 's': s_0, 'ds_dt': ds_dt}) 
# t_0_kep is not provided hence defaults to t_0

model_static = MulensModel.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho, 'q': q, 'alpha': alpha_0, 
     's': s_0})

# We can get model exactly the same as model_orb this way:
orb_parameters = {key: value for key, value in model_static.parameters.items()}
orb_parameters['ds_dt'] = 0.
orb_parameters['dalpha_dt'] = 0.
model_orb_2 = MulensModel.Model(Parameters=orb_parameters)

dt = 36.525 # This is in days.

################################################################
# Get the values of parameters in both models:
print(model_orb.s)
print(model_orb.s_0) # == model_orb.s
print(model_orb.get_s(t_0)) # Prints the same as previous one.
print(model_orb.get_s(t_0+dt)) # should return 1.55

print(model_static.s)
print(model_static.get_s(t_0)) # == model_static.s
print(model_static.get_s(t_0+dt)) # == model_static.s

print(model_orb.get_alpha(t_0-dt)) # should return 7.345 u.deg
# In analogy to s, similar methods for alpha will work.
################################################################

# Print projected orbital velocity
print(model_orb.gamma_parallel) # should return 0.3333333 1/u.year
print(model_orb.gamma_perp) # should return -50 1/u.year or -0.87266 u.rad/u.year 
# (the minus sign comes from the definition in Skowron et al. 2011)
print(model_orb.gamma) # should return 0.9346 1/u.year

# Make a nice plot:
plt.figure()
model_orb.plot_caustics(epoch=t_0)
model_orb.plot_caustics(epoch=t_0+dt, c='g') # second caustics are green
model_orb.plot_trajectory()
t = 'This plot shows a nice curved trajectory and caustics for 2 different epochs'
plt.title(t)
plt.show()

