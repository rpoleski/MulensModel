"""
Create Figure 1 of Poleski and Yee (2019)
"Modeling microlensing events with MulensModel"
Astronomy and Computing 26, 35
https://ui.adsabs.harvard.edu/abs/2019A&C....26...35P/abstract
https://arxiv.org/abs/1803.01003

Example point-source magnification curves.

"""
from matplotlib import pyplot

from MulensModel import Model


# Define model parameters.
t_0 = 2460000
params_ps = {'t_0': t_0, 'u_0': 0.0008, 't_E': 30.}
t_star = 0.051  # Day is default unit for t_E and t_star.
gamma = 0.4  # This is limb darkening coefficient.
params_fs = {**params_ps, 't_star': t_star}

# Set models settings for:
model_ps = Model(params_ps)  # point source,
model_fs = Model(params_fs)  # finite source,
model_fs_ld = Model(params_fs)  # and finite source with limb darkening.

# Specify which finite source methods are used and when:
t_1 = t_0 - 3.5 * t_star
t_2 = t_0 + 3.5 * t_star
model_fs.set_magnification_methods([t_1, 'finite_source_uniform_Gould94', t_2])
model_fs_ld.set_magnification_methods([t_1, 'finite_source_LD_Yoo04', t_2])

# Plot the magnification curves.
plot_kwargs = {'t_start': t_0-5.5*t_star, 't_stop': t_0+5.5*t_star,
               'subtract_2460000': True, 'lw': 2.}
model_ps.plot_magnification(label='point source', **plot_kwargs)
model_fs.plot_magnification(label='finite source', **plot_kwargs)
model_fs_ld.plot_magnification(
    gamma=gamma, label='finite source LD', **plot_kwargs)

pyplot.legend(loc='best')
pyplot.savefig('figure_1.png')
