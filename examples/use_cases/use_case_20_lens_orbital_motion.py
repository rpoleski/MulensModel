from astropy import units as u
import matplotlib.pyplot as plt

import MulensModel

"""
This use case needs a description.
"""

# point lens parameters:
t_0 = 2457123.456
u_0 = 0.0345
t_E = 30.00 * u.day
rho = 0.00345

# binary lens parameters:
q = 0.123
alpha = 12.345 * u.deg
dalpha_dt = 50. u.deg / u.year
s = 1.5
ds_dt = 0.5 / u.year

#Generate a model.
model = MulensModel.Model()
model.set_parameters(t_0=t_0, u_0=u_0, t_E=t_E, rho=rho, q=q, 
                    alpha=alpha, dalpha_dt=dalpha_dt, 
                    s=s, ds_dt=ds_dt, t_orb=t_0)
# QUESTION - should we have alpha_0 and s_0 instead of alpha and s?

# Get the values of parameters:
try:
    print(model.s) # Yes, this is not allowed. Same for alpha.
except: #this exception should be thrown from inside model.s. It should not be part of the use case. JCY
    print("You cannot call Model.s because lens orbital was specified " +
            "try Model.s_0 or Model.s_for_t(t))")
#OR
print(model.s) # should return input s value ( i.e. s_0 = s(t_orb) )
#RP and JCY disagree on the implementation

dt = 36.525 # This is in days.
print(model.s_for_epoch(t_0+dt)) # should return 1.55
print(model.alpha_for_epoch(epoch=t_0-dt)) # should return 7.345 u.deg

# Print projected orbital velocity
print(model.gamma_parallel) # should return 0.3333333 1/u.year
print(model.gamma_perp) # should return -50 1/u.year or -0.87266 u.rad/u.year 
# (the minus sign comes from the definition in Skowron et al. 2011)
print(model.gamma) # should return 0.9346 1/u.year

# Make a nice plot
plt.figure()
model.plot_caustics()
model.plot_caustics(epoch=t_0+dt)
model.plot_trajectory()
plt.title('This plot shows a nice curved trajectory')
plt.show()

