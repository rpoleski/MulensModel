"""
use_case_34_fix_blending.py

Fix the blending for one observatory to a non-zero value.

This example is modeled after real-time fitting procedures, so the
input datasets are truncated.

For this event, the catalog star in KMT is I > 20. Thus, a fake blending flux
has been added so the baseline is I=20. Thus, for a fixed blending fit, it is
appropriate to fix the blending at the value added to the baseline star.

Adapted from use_case_24_chi2_gradient.py
New functionality marked by # *** NEW ***
"""
import os
import MulensModel
import scipy.optimize as op


class Minimizer(object):
    """
    An object to link an Event to the functions necessary to minimize chi2.
    """

    def __init__(self, event, parameters_to_fit):
        self.event = event
        self.parameters_to_fit = parameters_to_fit

    def set_parameters(self, theta):
        """for given event set attributes from parameters_to_fit (list of str)
        to values from theta list"""
        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, val, theta[key])

    def chi2_fun(self, theta):
        """for a given set of parameters (theta), return the chi2"""
        self.set_parameters(theta)
        return self.event.get_chi2()

    def chi2_gradient(self, theta):
        """
        for a given set of parameters (theta), return the gradient of chi^2
        """
        self.set_parameters(theta)  # might be redundant, but probably safer
        return self.event.get_chi2_gradient(self.parameters_to_fit)


# *** NEW ***
datasets = []
file_names = ['KMTA12_I.pysis', 'KMTC12_I.pysis', 'KMTS12_I.pysis',
              'KMTA14_I.pysis', 'KMTC14_I.pysis', 'KMTS14_I.pysis']
data_ref = 1
dir_ = os.path.join(MulensModel.DATA_PATH, "photometry_files", "KB180003")
for i, file_name in enumerate(file_names):
    file_ = os.path.join(dir_, file_name)
    datasets.append(
        MulensModel.MulensData(
            file_name=file_, add_2450000=True, usecols=[0, 3, 4],
            phot_fmt='mag'))

# Calculate the amount of blending flux that was added.
Icat = 21.89
fcat = MulensModel.Utils.get_flux_from_mag(Icat)
fblend = MulensModel.Utils.get_flux_from_mag(20.) - fcat
# *** END NEW ***

# Fit t_0, t_E for several, fixed values of u_0 with
# fixed blending = added flux
# e.g. to assess whether or not the event could be high-magnification
parameters_to_fit = ["t_0", "t_E"]
for u_0 in [0.0, 0.01, 0.1]:
    t_0 = 2457215.
    t_E = 20.0
    model = MulensModel.Model(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    event = MulensModel.Event(datasets=datasets, model=model)
    # *** NEW ***
    event.fix_blend_flux[datasets[data_ref]] = fblend  # Fix the blending =
    # the known added value
    # *** END NEW ***

    initial_guess = [t_0, t_E]
    minimizer = Minimizer(event, parameters_to_fit)
    result = op.minimize(
        minimizer.chi2_fun, x0=initial_guess, method='Newton-CG',
        jac=minimizer.chi2_gradient, tol=1e-3)

    chi2 = minimizer.chi2_fun(result.x)
    print(chi2, u_0, result.x)
    print(minimizer.event.fluxes)  # *** NEW ***
