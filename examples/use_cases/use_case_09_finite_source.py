import MulensModel

#Delete? Redundant with UC18/19?

raise NotImplementedError('Model does not have a property finite_source')

model = MulensModel.Model(
    {'t_0': 7605.5, 'u_0': 0.001, 't_E': 1.0, 'rho': 0.001})
model.finite_source = True

model_2 = MulensModel.Model(
    {'t_0': 7605.5, 'u_0': 0.001, 't_E': 1.0, 'rho': 0.001})
model_2.finite_source(
    method='Stokes', time_range=[7605., 7606.], check_range=True, 
    tolerance=0.0001) # time_range and check_range are optional


