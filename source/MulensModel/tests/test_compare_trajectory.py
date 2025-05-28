import MulensModel as mm
import numpy as np
import unittest
import pytest

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

N = 5
times = []

for i in range(0,N,1):
    times.append(t_0 - 3 * t_E + i * (6 * t_E / (N - 1)))


parameters_classic = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho})

parameters_extra = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt})

def unit_change_for_VBB(ds_dt_MM, da_dt_MM, dsz_dt_MM):
    s_z = parameters_extra.s_z
    s_prim = np.sqrt(s**2 + s_z**2)
    
    ds_dt_VBB = ds_dt_MM * 1/s_prim * 1/365
    da_dt_VBB = da_dt_MM * np.pi/180 * 1/365
    dsz_dt_VBB = dsz_dt_MM * 1/s_prim * 1/365
    return ds_dt_VBB, da_dt_VBB, dsz_dt_VBB

#print(unit_change_for_VBB(ds_dt, dalpha_dt, ds_z_dt))

def test_init_parameters():
    """
    compare trajectory
    """
    trajectory = mm.Trajectory(parameters=parameters_extra, times=times)
    
    x_VBB = [3.026951, 1.597578, 0.166892, -1.265083, -2.698321]
    y_VBB = [0.381569, -0.075952, -0.529313, -0.978565, -1.423758]

    np.testing.assert_almost_equal(trajectory.x, x_VBB)
    np.testing.assert_almost_equal(trajectory.y, y_VBB)

