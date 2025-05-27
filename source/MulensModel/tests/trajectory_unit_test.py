import MulensModel as mm
import numpy as np
s = 1.2
q = 0.123
u_0 = 0.555
alpha = 17.5
t_0 = 2455500
t_E = 100
rho = 0.01


ds_dt = 0.013
dalpha_dt = -0.2
ds_z_dt = 0.05


np = 101
times = []

for i in range(0,np,1):
    times.append(t_0 - 3 * t_E + i * (6 * t_E / (np - 1)))


parameters_classic = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho})

parameters_extra = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt})

trajectory = mm.Trajectory(parameters=parameters_extra, times=times)
x = trajectory.x
y = trajectory.y


for i in range(0,np,1):
    print("t: %f  x: %f y: %f\n" %(times[i], x[i], y[i]));

