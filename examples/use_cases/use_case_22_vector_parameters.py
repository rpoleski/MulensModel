import MulensModel


# Show a new way of accessing model parameters.
t_0 = 2456789.0123
u_0 = 0.1
t_E = 15.
rho = 0.001 
pi_E_N = 0.1
pi_E_E = 0.2

parameters = MulensModel.ModelParameters(t_0=t_0, u_0=u_0, t_E=t_E, 
                pi_E_N=pi_E_N, pi_E_E=pi_E_E)
print(parameters.vector)
# returns np.array([t_0, u_0, t_E, pi_E_N, pi_E_E])

parameters.rho = rho
print(parameters.vector)
# returns np.array([t_0, u_0, t_E, rho, pi_E_N, pi_E_E])

