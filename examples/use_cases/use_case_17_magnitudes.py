"""
Setting source and blend fluxes for model plotting.
"""
import matplotlib.pyplot as plt

import MulensModel

# Define a Model
t_0 = 2456791.
u_0 = 0.2
t_E = 12.4
model_1 = MulensModel.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Plot the Model 3 ways.
plt.figure()
plt.title('Model Lightcurve defined using fluxes')
model_1.plot_lc(f_source=0.2, f_blend=0.4)

plt.figure()
plt.title('Model Lightcurve defined using magnitudes')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(21.2),
                f_blend=MulensModel.Utils.get_flux_from_mag(18.2))

plt.figure()
plt.title('Model Lightcurve in 2 bands')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(21.2),
                f_blend=MulensModel.Utils.get_flux_from_mag(18.2), label='I')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(20.3),
                f_blend=MulensModel.Utils.get_flux_from_mag(17.4), label='V')
plt.legend(loc='best')

plt.show()
