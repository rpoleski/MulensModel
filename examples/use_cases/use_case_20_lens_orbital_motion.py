"""
Use case presenting binary lens orbital motion models. Also comparison
with static binary model is provided.
"""
from astropy import units as u
import matplotlib.pyplot as plt

from MulensModel import Model


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
params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho, 'q': q, 
        'alpha': alpha_0, 's': s_0}

#Generate models:
model_orb = Model({**params, 'dalpha_dt': dalpha_dt, 'ds_dt': ds_dt})
# t_0_kep is not provided hence defaults to t_0
model_static = Model(params)

# We can get model exactly the same as model_orb this way:
orb_parameters = model_static.parameters.as_dict().copy()
orb_parameters.update({'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt})
model_orb_2 = MulensModel.Model(Parameters=orb_parameters)

dt = 36.525 # This is in days.

# Get the values of parameters in both models:
print(model_orb.parameters.s)
print(model_orb.parameters.get_s(t_0)) # Prints the same as previous one.
print(model_orb.parameters.get_s(t_0+dt)) # should return 1.55

print(model_static.parameters.s)
print(model_static.parameters.get_s(t_0)) # == model_static.s
print(model_static.parameters.get_s(t_0+dt)) # == model_static.s

print(model_orb.parameters.get_alpha(t_0-dt)) # should return 7.345 u.deg
# In analogy to s, similar methods for alpha will work.

# Make sure that you know what kind of model you deal with:
print(model_static.is_static()) # Returns True.
print(model_orb.is_static()) # Returns False.
print(model_orb_2.is_static()) # Returns False.

# Print projected orbital velocity
print(model_orb.parameters.gamma_parallel) # Should return 0.3333333 1/u.year
print(model_orb.parameters.gamma_perp) # Should return -50 1/u.year or 
# -0.87266 u.rad/u.year (the minus sign comes from the definition in 
# Skowron et al. 2011).
print(model_orb.parameters.gamma) # Should return 0.9346 1/u.year

# Make a nice plot:
plt.figure()
model_orb.plot_caustics(epoch=t_0)
model_orb.plot_caustics(epoch=t_0+dt, c='g') # second caustics are green
model_orb.plot_trajectory()
plt.title('This plot shows a nice curved trajectory and caustics for 2 ' + 
        'different epochs')
plt.show()

