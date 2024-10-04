"""
Use case presenting binary lens orbital motion models. Also comparison
with static binary model is provided.
"""
import matplotlib.pyplot as plt

import MulensModel as mm


# point lens parameters:
t_0 = 2457123.456
u_0 = 0.0345
t_E = 30.00
rho = 0.00345

# binary lens parameters:
q = 0.123
alpha_0 = 12.345
dalpha_dt = 50.
s_0 = 1.5
ds_dt = 0.5
params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho, 'q': q,
          'alpha': alpha_0, 's': s_0}

# Generate models:
model_static = mm.Model(params)
params['dalpha_dt'] = dalpha_dt
params['ds_dt'] = ds_dt
model_orb = mm.Model(params)
# t_0_kep is not provided hence defaults to t_0

# We can get model exactly the same as model_orb this way:
orb_parameters = model_static.parameters.as_dict().copy()
orb_parameters.update({'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt})
model_orb_2 = mm.Model(parameters=orb_parameters)

dt = 36.525  # This is in days.

# Get the values of parameters in both models:
print('test orbital motion s')
print('{0} == {1}?'.format(model_orb.parameters.s, s_0))
print('{0} == {1}?'.format(model_orb.parameters.get_s(t_0), s_0))
print('{0} == {1}?'.format(model_orb.parameters.get_s(t_0+dt), 1.55))

print('test static s')
print('{0} == {1}?'.format(model_static.parameters.s, s_0))
print('{0} == {1}?'.format(model_static.parameters.get_s(t_0), s_0))
print('{0} == {1}?'.format(model_static.parameters.get_s(t_0+dt), s_0))

print('test orbital motion alpha')
print(
    '{0} == {1}?'.format(
        model_orb.parameters.get_alpha(t_0-dt), 7.345))
# In analogy to s, similar methods for alpha will work.

# Make sure that you know what kind of model you deal with:
print('What kind of model?')
print('Static: {0} == {1}?'.format(model_static.is_static(), True))
print('Orb Mot: {0} == {1}?'.format(model_orb.is_static(), False))
print('Orb Mot: {0} == {1}?'.format(model_orb_2.is_static(), False))

# Print projected orbital velocity
print('Check projected velocity (gamma)')
print(
    '{0} == {1}?'.format(
        model_orb.parameters.gamma_parallel, 0.3333333))
print(
    '{0} == {1}?'.format(
        model_orb.parameters.gamma_perp,  -0.87266))
# or # -0.87266 u.rad/u.year (the minus sign comes from the definition in
# Skowron et al. 2011).
print('{0} == {1}?'.format(model_orb.parameters.gamma, 0.9346))

# Make a nice plot:
plt.figure()
model_orb.plot_caustics(epoch=t_0)
model_orb.plot_caustics(epoch=t_0+dt, c='g', lw=0)  # second caustics are green
model_orb.plot_trajectory()
plt.title('This plot shows a nice curved trajectory and caustics for 2 ' +
          'different epochs')
plt.show()
