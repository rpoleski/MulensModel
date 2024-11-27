"""
Plots 2 light curves that are intended to be an inspiration for
a MulensModel logo.
"""
import MulensModel as mm
import matplotlib.pyplot as plt

binary_model = mm.Model(
    {'t_0': 1.3, 'u_0': 0.11, 't_E': 1.5, 's': 0.5, 'q': 0.3, 'alpha': 142.})

binary_source_model = mm.Model(
    {'t_0_1': 0.5, 'u_0_1': 0.08, 't_E': 0.5,
     't_0_2': 0.8, 'u_0_2': 0.11})

plt.figure()
binary_source_model.plot_lc(
    source_flux=[1.4, 1.6], blend_flux=18.0, color='blue', lw=10)
binary_model.plot_lc(
    source_flux=1.0, blend_flux=20.0, color='magenta', lw=10)
plt.xlim(0., 2.)
plt.show()
