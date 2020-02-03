import numpy as np
from numpy.testing import assert_almost_equal as almost
import unittest

import MulensModel as mm


def generate_model():
    """
    returns a model, time array, and magnification
    """

    # Create a PSPL model
    t_0 = 3583.
    u_0 = 0.3
    t_E = 12.

    t = np.linspace(t_0 - 3. * t_E, t_0 + 3. * t_E, 1000)
    pspl = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    A = pspl.magnification(t)

    return (pspl, t, A)


def generate_binary_model():
    """
    returns a binary source model, time array, and the magnification of
    both sources
    """

    # retrieve model 1
    (model_1, t, A_1) = generate_model()
    t_0_1 = model_1.parameters.t_0
    u_0_1 = model_1.parameters.u_0

    # create second model
    t_0_2 = 3570.
    u_0_2 = 0.25
    t_E = model_1.parameters.t_E

    model_2 = mm.Model({'t_0': t_0_2, 'u_0': u_0_2, 't_E': t_E})

    A_2 = model_2.magnification(t)

    # create separate binary model
    params = {'t_0_1': t_0_1, 'u_0_1': u_0_1, 't_0_2': t_0_2, 'u_0_2': u_0_2,
              't_E': t_E}
    binary_model = mm.Model(params)

    return (binary_model, t, A_1, A_2)


def generate_dataset(f_mod, t):
    """
    pass in f_mod and t, returns a MulensData
    """

    # error in measurement
    err = f_mod * 0.01

    my_dataset = mm.MulensData(data_list=[t, f_mod, err], phot_fmt='flux')

    return my_dataset


def execute_test_blend_fixed(f_b):
    # test for when source flux is to be determined, but blend flux is a
    # fixed value

    pspl, t, A = generate_model()

    # secret source flux
    f_s = 1.0
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=f_b,
        fix_source_flux=False)
    my_fit.fit_fluxes()

    almost(my_fit.source_flux, f_s)


class BinarySourceTest():

    def __init__(self):
        self.f_s_1 = 1
        self.f_s_2 = 1.2
        self.f_b = 0.5

        self.setup_dataset()

    def setup_dataset(self):
        self.model, self.t, self.A_1, self.A_2 = generate_binary_model()
        f_mod = self.f_s_1 * self.A_1 + self.f_s_2 * self.A_2 + self.f_b
        self.dataset = generate_dataset(f_mod, self.t)

    def run_test(
            self, fix_blend_flux=False, fix_source_flux=False,
            fix_q_flux=False):

        my_fit = mm.FitData(
            model=self.model, dataset=self.dataset,
            fix_blend_flux=fix_blend_flux,
            fix_source_flux=fix_source_flux, fix_q_flux=fix_q_flux)
        my_fit.fit_fluxes()

        almost(my_fit.blend_flux, self.f_b)
        almost(my_fit.source_fluxes[0], self.f_s_1)
        almost(my_fit.source_fluxes[1], self.f_s_2)


def execute_test_binary_source(q_flux=False):
    # test for when blend flux and source flux are to be determined for binary
    # sources with q-flux

    test = BinarySourceTest()
    if q_flux:
        fix_q_flux = test.f_s_2 / test.f_s_1
    else:
        fix_q_flux = False

    test.run_test(fix_q_flux=fix_q_flux)


# *** Actual tests below ***
def test_default():
    """
    test for when blend flux and source flux are to be determined
    """
    pspl, t, A = generate_model()

    # secrets
    f_s = 1.0
    f_b = 0.5
    # generate f_mod
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=False)
    my_fit.fit_fluxes()

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_flux, f_s)


def test_blend_zero():
    """
    test for when source flux is to be determined, but blend flux is zero
    """
    execute_test_blend_fixed(f_b=0.)


def test_blend_fixed():
    """
    test for when source flux is to be determined, but blend flux is
    zero
    """
    execute_test_blend_fixed(f_b=0.5)


def test_source_fixed():
    """
    test for when blend flux is to be determined, but source flux is a fixed
    value
    """

    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1.0
    f_b = 0.5
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=f_s)
    my_fit.fit_fluxes()

    almost(my_fit.blend_flux, f_b)


def test_binary_source():
    """Test a binary source model with all free parameters."""
    execute_test_binary_source(q_flux=False)


def test_binary_source_fixed():
    """
    Test the three cases for fixing each of the three components of a binary
    source model
    """
    test = BinarySourceTest()
    test.run_test(fix_source_flux=[1.0, False])
    test.run_test(fix_source_flux=[False, 1.2])
    test.run_test(fix_source_flux=[1.0, 1.2])
    test.run_test(fix_blend_flux=0.5)


class TestFitData(unittest.TestCase):
    def test_init_1(self):
        with self.assertRaises(ValueError):
            test = BinarySourceTest()
            test.run_test(fix_source_flux=1.0)


def test_binary_qflux():
    """
    test for when blend flux and source flux are to be determined for binary
    sources with q-flux
    """

    execute_test_binary_source(q_flux=True)

def test_fit_fluxes():
    """
    test that when the model is updated, and fit fluxes is re-run, the fluxes
    actually change.
    """

    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1.0
    f_b = 0.5
    f_mod = f_s * A + f_b

    my_dataset = generate_dataset(f_mod, t)
    my_fit = mm.FitData(
        model=pspl, dataset=my_dataset, fix_blend_flux=False,
        fix_source_flux=False)
    my_fit.update()
    f_s_1 = my_fit.source_flux
    chi2_1 = my_fit.chi2

    t_E_2 = pspl.parameters.t_E / (f_s + f_b)
    u_0_2 = pspl.parameters.u_0 / (f_s + f_b)
    new_model = mm.Model(
        {'t_0': pspl.parameters.t_0, 'u_0': u_0_2, 't_E': t_E_2})
    my_fit.model = new_model
    my_fit.fix_blend_flux = 0.
    my_fit.fit_fluxes()

    assert(f_s_1 != my_fit.source_flux)
    assert(chi2_1 == my_fit.chi2)

    my_fit.update()
    assert(chi2_1 != my_fit.chi2)

def test_chi2_per_point():
    """Test that the chi2 shape is correct for multiple sources, i.e. = number
    of epochs, rather than epochs * sources."""
    test_object = BinarySourceTest()
    my_fit = mm.FitData(model=test_object.model, dataset=test_object.dataset)
    my_fit.update()

    assert(my_fit.chi2_per_point.shape == (test_object.dataset.n_epochs,))

