"""
Note that this use case is not consistent with any other use_cases. It
is meant to work out how this should be done. This definition should
be retroactively applied to other use cases.
"""
import matplotlib.pyplot as pl

import MulensModel

t_0 = 6791.
u_0 = 0.2
t_E = 12.4

model_1 = MulensModel.Model(t_0=t_0, u_0=u_0, t_E=t_E, f_s=0.2, f_b=0.4)
pl.plot(model_1.time, model_1.magnitude)

model_2 = MulensModel.Model(t_0=t_0, u_0=u_0, t_E=t_E)
model_2.source.I_mag = 21.2
model_2.blend.I_mag = 18.2
pl.plot(model_2.time, model_2.magnitude)

model_3 = MulensModel.Model(t_0=t_0, u_0=u_0, t_E=t_E)
model_3.source.magnitude(21.2, bandpass='I')
model_3.blend.magnitude(18.2, bandpass='I')
model_3.source.magnitude(20.3, bandpass='V')
model_3.blend.magnitude(17.4, bandpass='V')
pl.plot(model_3.time, model_3.magnitude) #defaults to I?
pl.plot(model_3.time, model_3.V_mag)
