import numpy as np

import MulensModel as mm


def test_magnification_type():
    """
    Check type of magnification returned for model with t_eff.
    At some point it was astropy quantity.
    """
    parameters = mm.ModelParameters({'t_0': 1., 't_eff': 0.2, 't_E': 3.})
    magnification_curve = mm.MagnificationCurve(2., parameters)
    assert type(magnification_curve.get_magnification()) == np.ndarray


def test_fspl_noLD():
    """
    check if FSPL magnification is calculate properly
    """
    t_0 = 2456789.012345
    t_E = 23.4567
    u_0 = 1e-4
    rho = 1e-3
    t_vec = np.array([-(rho**2-u_0**2)**0.5, 0., ((0.5*rho)**2-u_0**2)**0.5])
    t_vec = t_vec * t_E + t_0

    params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})

    mag_curve = mm.MagnificationCurve(times=t_vec, parameters=params)
    methods = [t_0-t_E, 'finite_source_uniform_Gould94', t_0+t_E]
    mag_curve.set_magnification_methods(methods, 'point_source')
    results = mag_curve.get_point_lens_magnification()

    u = np.array([rho, u_0, 0.5*rho])
    pspl = (u**2 + 2.) / np.sqrt(u**2 * (u**2 + 4.))
    expected = np.array([1.27323965, 0.19949906, 0.93421546])
# These values were calculated by Andy Gould (file b0b1.dat).
    expected *= pspl

    np.testing.assert_almost_equal(expected, results, decimal=4)


def test_fspl():
    """
    check if FSPL magnification is calculate properly
    """
    t_0 = 2456789.012345
    t_E = 23.4567
    u_0 = 1e-4
    rho = 1e-3
    gamma = 0.56789
    t_vec = np.array([-(rho**2-u_0**2)**0.5, 0., ((0.5*rho)**2-u_0**2)**0.5])
    t_vec = t_vec * t_E + t_0

    params = mm.ModelParameters(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})

    mag_curve = mm.MagnificationCurve(
        times=t_vec, parameters=params, gamma=gamma)
    methods = [t_0-t_E, 'finite_source_LD_Yoo04', t_0+t_E]
    mag_curve.set_magnification_methods(methods, 'point_source')
    results = mag_curve.get_point_lens_magnification()

    u = np.array([rho, u_0, 0.5*rho])
    pspl = (u**2 + 2.) / np.sqrt(u**2 * (u**2 + 4.))
    expected = np.array([1.27323965-gamma*0.09489869,
                         0.19949906-gamma*-0.03492121,
                         0.93421546-gamma*-0.09655794])
# These values were calculated by Andy Gould (file b0b1.dat).
    expected *= pspl

    np.testing.assert_almost_equal(expected/results, 1., decimal=4)


def test_Lee09_and_WittMao94():
    """
    test Lee et al. 2009 and Witt & Mao 1994 finite source calculation
    """
    t_vec = np.array([3.5, 2., 1., 0.5, 0.])

# The values below were calculated using code developed by P. Mroz.
    expected_0 = np.array([1.01084060513, 1.06962639343, 1.42451408166,
                           2.02334097551, 2.13919086656])
    expected_1 = np.array([1.01110609638, 1.07461016241, 1.57232954942,
                           2.21990790526, 2.39458814753])
    expected_2 = np.array([1.0110829794, 1.07404148634, 1.55620547462,
                           2.24809136704, 2.44503143812])
# The last values are for 2-parameter LD with same settings and lambda=0.3.
# Correction is:
#  -lambda*(1-1.25*sqrt(costh))
# and for 1-parameter LD we used:
#  1-gamma*(1-1.5*costh)

    # Test uniform source first.
    params_0 = mm.ModelParameters(
        {'t_0': 0., 'u_0': 0.5, 't_E': 1., 'rho': 1.})
    mag_curve_0 = mm.MagnificationCurve(times=t_vec, parameters=params_0)
    methods_0 = [-5., 'finite_source_uniform_Lee09', 5.]
    mag_curve_0.set_magnification_methods(methods_0, 'point_source')
    results_0 = mag_curve_0.get_point_lens_magnification()
    np.testing.assert_almost_equal(expected_0, results_0, decimal=4)

    # Then test 1-parameter limb-darkening.
    params_1 = mm.ModelParameters(
        {'t_0': 0., 'u_0': 0.1, 't_E': 1., 'rho': 1.})
    mag_curve_1 = mm.MagnificationCurve(times=t_vec, parameters=params_1,
                                        gamma=0.5)
    methods_1 = [-5., 'finite_source_LD_Lee09', 5.]
    mag_curve_1.set_magnification_methods(methods_1, 'point_source')
    results_1 = mag_curve_1.get_point_lens_magnification()
    np.testing.assert_almost_equal(expected_1, results_1, decimal=3)

    # Tests for Witt & Mao 1994 start here
    methods_2 = [-5., 'finite_source_uniform_WittMao94', 5.]
    mag_curve_0.set_magnification_methods(methods_2, 'point_source')
    results_2 = mag_curve_0.get_point_lens_magnification()
    np.testing.assert_almost_equal(expected_0, results_2, decimal=4)

    methods_3 = [-5., 'finite_source_LD_WittMao94', 5.]
    mag_curve_1.set_magnification_methods(methods_3, 'point_source')
    results_3 = mag_curve_1.get_point_lens_magnification()
    np.testing.assert_almost_equal(expected_1, results_3, decimal=3)


def test_PSPL_for_binary():
    """
    test PSPL model used in a model that is defined as binary
    """
    t_0 = 1000.
    t_E = 20.
    u_0 = 1.
    t_vec = np.array([10., 100.]) * t_E + t_0
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 's': 1.2, 'q': 0.1, 'alpha': 0.})
    mag_curve = mm.MagnificationCurve(times=t_vec, parameters=params)
    mag_curve.set_magnification_methods(None, 'point_source_point_lens')
    u2 = u_0**2 + ((t_vec - t_0) / t_E)**2
    pspl = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    np.testing.assert_almost_equal(pspl, mag_curve.get_magnification())
