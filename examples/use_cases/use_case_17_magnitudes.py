"""
Note that this use case is not consistent with any other use_cases. It
is meant to work out how this should be done. This definition should
be retroactively applied to other use cases.
"""
import matplotlib.pyplot as pl

import MulensModel

#Define a Model
t_0 = 2456791.
u_0 = 0.2
t_E = 12.4
model_1 = MulensModel.Model(t_0=t_0, u_0=u_0, t_E=t_E)

#Plot the Model 3 ways.
pl.figure()
pl.title('Model Lightcurve defined using fluxes')
model_1.plot_lc(f_source=0.2, f_blend=0.4)

pl.figure()
pl.title('Model Lightcurve defined using magnitudes')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(21.2),
                f_blend=MulensModel.Utils.get_flux_from_mag(18.2))

pl.figure()
pl.title('Model Lightcurve in 2 bands')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(21.2),
                f_blend=MulensModel.Utils.get_flux_from_mag(18.2),label='I')
model_1.plot_lc(f_source=MulensModel.Utils.get_flux_from_mag(20.3),
                f_blend=MulensModel.Utils.get_flux_from_mag(17.4),label='V')
pl.legend(loc='best')

pl.show()
