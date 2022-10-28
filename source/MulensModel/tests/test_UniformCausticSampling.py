import numpy as np

import MulensModel as mm


def test_UniformCausticSampling_simple():
    """
    Tests only very basic usage of UniformCausticSampling.
    """
    sampling = mm.UniformCausticSampling(s=1.5, q=0.9, n_points=10)
    assert sampling.s == 1.5
    assert sampling.q == 0.9
    assert sampling.n_caustics == 1


def test_UniformCausticSampling():
    """
    Test if the main functions of work correctly, namely:
    - get_standard_parameters(),
    - get_x_in_x_out(),
    - check_valid_trajectory(),
    - caustic_point(),
    - n_caustics,
    - which_caustic().
    Most tests are run for each of the three topologies.
    """
    s_ = [1.1, 0.5, 2]
    q_ = [0.1, 0.001, 0.01]
    n_points = 1000
    # Expected results:
    n_c = [1, 3, 2]
    p_ = [0.01401270+0.22985537j, -1.49730846-0.11183662j,
          1.49051495+0.03108463j]
    checks = [False, False, True]

    for (s, q, nc, p, c) in zip(s_, q_, n_c, p_, checks):
        sampling = mm.UniformCausticSampling(s=s, q=q, n_points=n_points)
        standard = sampling.get_standard_parameters(0.2, 0.499, 0., 10.)
        x_list = sampling.get_x_in_x_out(standard['u_0'], standard['alpha'])
        assert np.min(np.abs(np.array(x_list) - 0.2)) < 1.e-6
        assert np.min(np.abs(np.array(x_list) - 0.499)) < 2.e-6
        assert sampling.check_valid_trajectory(0.2, 0.5) == c
        np.testing.assert_almost_equal(sampling.caustic_point(0.3), p)
        assert sampling.n_caustics == nc
        assert sampling.which_caustic(0.95) == sampling.n_caustics
    np.testing.assert_almost_equal(standard['t_0'], 216.3783662)
    np.testing.assert_almost_equal(standard['u_0'], 0.0076099741)
    np.testing.assert_almost_equal(standard['t_E'], 142.7154375)
    np.testing.assert_almost_equal(standard['alpha'], 180.54305134)

    # For the last (s, q) we also check alpha = 90 or 270:
    x_1 = 0.1885475223517683
    x_2 = 0.9504818975429649
    u_0 = s - 1. / s + s * q / (1. + q)
    standard = sampling.get_standard_parameters(x_1, x_2, -10., 0.)
    np.testing.assert_almost_equal(standard['alpha'], 270.)
    np.testing.assert_almost_equal(standard['u_0'], u_0)
    standard = sampling.get_standard_parameters(x_2, x_1, -10., 0.)
    np.testing.assert_almost_equal(standard['alpha'], 90.)
    np.testing.assert_almost_equal(standard['u_0'], -u_0)


def test_get_uniform_sampling():
    """
    Basic test for get_uniform_sampling() - only if the shape of output is
    correct
    """
    n = 50

    sampling = mm.UniformCausticSampling(s=1.1, q=0.5, n_points=100)
    (x_caustic_in, x_caustic_out) = sampling.get_uniform_sampling(n_points=n)
    assert x_caustic_in.shape == (n,)
    assert x_caustic_out.shape == (n,)
