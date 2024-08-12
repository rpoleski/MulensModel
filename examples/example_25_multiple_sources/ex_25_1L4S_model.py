"""
Plot a 1L4S  model. Also show the __repr__ string for that model.
"""
import MulensModel as mm
import matplotlib.pyplot as plt

t_E = 10.
source_1 = {'t_0': 0., 'u_0': 0.1, 'source_flux': 1.}
source_2 = {'t_0': -5., 'u_0': 0.3, 'source_flux': 2.}
source_3 = {'t_0': 1., 'u_0': 0.001, 'rho': 0.002, 'source_flux': 0.01}
source_4 = {'t_0': 1.4, 'u_0': 0.0003, 'rho': 0.0005, 'source_flux': 0.005}

model_params = {'t_E': t_E}
source_fluxes = []
models = []
for i, source in zip(range(1, 5), [source_1, source_2, source_3, source_4]):
    single_source_model_params = {'t_E': t_E}
    for key, value in source.items():
        if key in ['t_0', 'u_0', 'rho']:
            model_params['{0}_{1}'.format(key, i)] = value
            single_source_model_params[key] = value
        elif key == 'source_flux':
            source_fluxes.append(value)

    models.append(mm.Model(single_source_model_params))

print('dict model_params:\n', model_params)

model = mm.Model(model_params)
print('model.__repr__:\n', model)

model.plot_lc(source_flux=source_fluxes, label='1L4S Model', color='black')
for i, single_model, source_flux in zip(range(1, 5), models, source_fluxes):
    single_model.plot_lc(
        source_flux=source_flux, label='Source {0}'.format(i))

plt.title('1L4S model + Individual Components')
plt.legend()
plt.minorticks_on()

plt.figure()
plt.gca().set_aspect('equal')
plt.title('Source Trajectories')
model.plot_trajectory()
plt.scatter(0, 0, marker='x', color='black', zorder=10)
plt.minorticks_on()
plt.xlim(-1.5, 0.5)
plt.ylim(-1., 1.)

plt.show()
