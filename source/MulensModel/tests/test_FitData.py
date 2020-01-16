import numpy as np
from numpy.testing import assert_almost_equal as almost

import MulensModel as mm


def generate_model():
    """
    returns a model, time array, and magnification
    """

    # Create a PSPL model
    t_0 = 3583.
    u_0 = 0.3
    t_E = 12.

    pspl = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

    t = np.linspace(3576, 3590, 1000)

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

    my_dataset = mm.MulensData(data_list = [t, f_mod, err], phot_fmt='flux')

    return my_dataset

def test_default():
    # test for when blend flux and source flux are to be determined
    
    pspl, t, A = generate_model()

    # secrets
    f_s = 1
    f_b = 0.5

    # generate f_mod
    f_mod = f_s * A + f_b 

    my_dataset = generate_dataset(f_mod, t)

    my_fit = mm.FitData(model = pspl, dataset = my_dataset,
        fix_blend_flux=False, fix_source_flux=False)
    
    my_fit.fit_fluxes()

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_flux, f_s)

def test_blend_zero():
    # test for when source flux is to be determined, but blend flux is zero
    test_blend_fixed(f_b = 0)

def test_blend_fixed(f_b = 0.5):
    # test for when source flux is to be determined, but blend flux is a fixed value

    pspl, t, A = generate_model()

    # secret source flux
    f_s = 1

    f_mod = f_s * A + f_b 

    my_dataset = generate_dataset(f_mod, t)

    my_fit = mm.FitData(model = pspl, dataset = my_dataset, f_blend = f_b, f_source = True)

    almost(my_fit.source_flux, f_s)

def test_source_fixed():
    # test for when blend flux is to be determined, but source flux is a fixed value

    pspl, t, A = generate_model()

    # secret blend flux, set source flux
    f_s = 1
    f_b = 0.5

    f_mod = f_s * A + f_b 

    my_dataset = generate_dataset(f_mod, t)

    my_fit = mm.FitData(model = pspl, dataset = my_dataset, f_blend = True, f_source = f_s)

    almost(my_fit.blend_flux, f_b)

def test_binary(q_flux):
    # test for when blend flux and source flux are to be determined for binary sources

    model, t, A_1, A_2 = generate_binary_model()

    # secrets
    f_s_1 = 1
    f_s_2 = 1.2
    f_b = 0.5

    f_mod = f_s_1 * A_1 + f_s_2 * A_2 + f_b

    my_dataset = generate_dataset(f_mod, t)

    my_fit = mm.FitData(model = model, dataset = my_dataset, fix_blend_flux=True, fix_source_flux=True)

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_fluxes[0], f_s_1)
    almost(my_fit.source_fluxes[1], f_s_2)


def test_binary_qflux():
    # test for when blend flux and source flux are to be determined for binary sources with q-flux

    model, t, A_1, A_2 = generate_binary_model()

    # secrets
    f_s_1 = 1
    f_s_2 = 1.2
    f_b = 0.5

    f_mod = f_s_1 * A_1 + f_s_2 * A_2 + f_b

    my_dataset = generate_dataset(f_mod, t)

    my_fit = mm.FitData(model = model, dataset = my_dataset, f_blend = True, f_source = True, q_flux = f_s_2/f_s_1)

    almost(my_fit.blend_flux, f_b)
    almost(my_fit.source_fluxes[0], f_s_1)
    almost(my_fit.source_fluxes[1], f_s_2)