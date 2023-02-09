"""
Plots various models for MB08310 using different implementations of
limb-darkening coefficients.

From `Janczak et al. 2010, ApJ 711, 731
<https://ui.adsabs.harvard.edu/abs/2010ApJ...711..731J/abstract>`_.

"""

import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm

# Define basic point lens model
t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
parameters = {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star}
(t_fs_min, t_fs_max) = (t_0 - 2. * t_star, t_0 + 2. * t_star)

# Uniform source model
uniform_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star},
    magnification_methods=[t_fs_min, 'finite_source_uniform_Gould94', t_fs_max],
    plot_properties={'label': 'uniform', 'color': 'black', 'linestyle': ':',
                     'zorder': 10})

# Yoo04 limb-darkening with gamma
# gamma_U = (gamma_I + gamma_V) / 2.
yoo_gamma_ld_coeffs = mm.LimbDarkeningCoeffs(gammas={'I': 0.4390, 'U': 0.526})
yoo_gamma_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star},
    magnification_methods=[t_fs_min, 'finite_source_LD_Yoo04', t_fs_max],
    limb_darkening_coeffs=yoo_gamma_ld_coeffs,
    plot_properties={'label': 'Yoo04'})

# Yoo04 limb-darkening with u
# From Claret 1998 for T=5800K, logg=5.0
# u_V = 0.704
yoo_u_ld_coeffs = mm.LimbDarkeningCoeffs(u={'I': 0.540, 'U': 0.625})
yoo_u_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star},
    magnification_methods=[t_fs_min, 'finite_source_LD_Yoo04', t_fs_max],
    limb_darkening_coeffs=yoo_u_ld_coeffs,
    plot_properties={'label': 'Yoo04-u'})

# 2-parameter limb-darkening with (gamma, lambda)
an_gamlam_ld_coeffs = mm.TwoParamLimbDarkeningCoeffs(
    gammas={'I': (0.077, 0.549), 'U': (0.166, 0.543)})
an_gamlam_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star},
    magnification_methods=[t_fs_min, 'finite_source_2paramLD_An02', t_fs_max],
    limb_darkening_coeffs=an_gamlam_ld_coeffs,
    plot_properties={'label': 'An02'})


# 2-parameter limb-darkening with (c, d)
an_cd_ld_coeffs = mm.TwoParamLimbDarkeningCoeffs(
    cds={'I': (0.099, 0.584), 'U': (0.204, 0.557)})
an_cd_model = mm.Model(
    {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star},
    magnification_methods=[t_fs_min, 'finite_source_2paramLD_An02', t_fs_max],
    limb_darkening_coeffs=an_cd_ld_coeffs,
    plot_properties={'label': 'An02-cd'})


# Time vectors for plotting (from Ex 05)
t_start = t_0 - 3.
t_stop = t_0 + 1.
n_star = 2.
t_star_start = t_0 - n_star * t_star
t_star_stop = t_0 + n_star * t_star
times = np.arange(t_start, t_star_start, 0.01)
times = np.concatenate((times, np.arange(t_star_start, t_star_stop, 0.0001)))
times = np.concatenate((times, np.arange(t_star_stop, t_stop, 0.01)))

# Make some plots
def common_plot_elements(title):
    plt.title(title)
    uniform_model.plot_magnification(times, subtract_2450000=True)
    plt.xlim(t_start - 2450000., t_stop - 2450000.)
    plt.legend(loc='upper left')

def plot_LD_model(model, bandpass=None, color=None):
    linestyle = {'I': '-', 'U': '-.'}
    model.plot_magnification(
        times, subtract_2450000=True,
        bandpass=bandpass, color=color, linestyle=linestyle[bandpass],
        label='{0} {1}'.format(model.plot_properties['label'], bandpass))

plt.figure(figsize=[8, 8])
plt.subplot(2, 2, 1)
common_plot_elements(r'Yoo04, $\gamma$')
plot_LD_model(yoo_gamma_model, bandpass='I', color='magenta')
plot_LD_model(yoo_gamma_model, bandpass='U', color='blue')

plt.subplot(2, 2, 2)
common_plot_elements(r'An02, $(\Gamma, \Lambda)$')
plot_LD_model(an_gamlam_model, bandpass='I', color='red')
plot_LD_model(an_gamlam_model, bandpass='U', color='cyan')

plt.subplot(2, 2, 3)
common_plot_elements(r'I: Yoo04, $u$ vs. An02, $(\Gamma, \Lambda)$')
plot_LD_model(yoo_gamma_model, bandpass='I', color='magenta')
plot_LD_model(an_gamlam_model, bandpass='I', color='red')

plt.subplot(2, 2, 4)
common_plot_elements(r'U: Yoo04, $u$ vs. An02, $(\Gamma, \Lambda)$')
plot_LD_model(yoo_gamma_model, bandpass='U', color='blue')
plot_LD_model(an_gamlam_model, bandpass='U', color='cyan')

plt.show()