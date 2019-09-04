import numpy as np

from MulensModel import UniformCausticSampling


def test_UniformCausticSampling_simple():
    """
    Tests only very basic usage of UniformCausticSampling.
    """
    sampling = UniformCausticSampling(s=1.5, q=0.9, n_points=10)
    assert sampling.s == 1.5
    assert sampling.q == 0.9
    assert sampling.n_caustics == 1


def test_UniformCausticSampling():
    """
    Test if the main functions of work correctly, namely:
    - get_standard_parameters(),
    - get_x_in_x_out(),
    - orientation_check(),
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
    p_ = [0.0140127+0.2298554j, -1.4973085-0.1118366j, 1.4905150+0.0310846j]
    checks = [False, False, True]

    for (s, q, nc, p, c) in zip(s_, q_, n_c, p_, checks):
        sampling = UniformCausticSampling(s=s, q=q, n_points=n_points)
        standard = sampling.get_standard_parameters(0.2, 0.5, 0., 10.)
        x_list = sampling.get_x_in_x_out(standard['u_0'], standard['alpha'])
        assert np.min(np.abs(np.array(x_list) - 0.2)) < 1.e-6
        assert np.min(np.abs(np.array(x_list) - 0.5)) < 1.e-6
        assert sampling.orientation_check(0.2, 0.5) == c
        np.testing.assert_almost_equal(sampling.caustic_point(0.3), p)
        assert sampling.n_caustics == nc
        assert sampling.which_caustic(0.95) == sampling.n_caustics
    np.testing.assert_almost_equal(standard['t_0'], 215.3545535)
    np.testing.assert_almost_equal(standard['u_0'], 0.0104627)
    np.testing.assert_almost_equal(standard['t_E'], 142.0417601)
    np.testing.assert_almost_equal(standard['alpha'], 180.6508586)


def test_get_uniform_sampling():
    """
    Basic test for get_uniform_sampling() - only if the shape of output is
    correct
    """
    n = 50

    sampling = UniformCausticSampling(s=1.1, q=0.5, n_points=100)
    (x_caustic_in, x_caustic_out) = sampling.get_uniform_sampling(n_points=n)
    assert x_caustic_in.shape == (n,)
    assert x_caustic_out.shape == (n,)
