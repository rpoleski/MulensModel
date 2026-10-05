import unittest

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
    np.testing.assert_almost_equal(standard['alpha'], 0.54305134)

    # For the last (s, q) we also check alpha = 90 or 270:
    x_1 = 0.1885475223517683
    x_2 = 0.9504818975429649
    u_0 = s - 1. / s + s * q / (1. + q)
    standard = sampling.get_standard_parameters(x_1, x_2, -10., 0.)
    np.testing.assert_almost_equal(standard['alpha'], 90.)
    np.testing.assert_almost_equal(standard['u_0'], u_0)
    standard = sampling.get_standard_parameters(x_2, x_1, -10., 0.)
    np.testing.assert_almost_equal(standard['alpha'], 270.)
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


class TestSelectNPoints(unittest.TestCase):
    """
    Tests for UniformCausticSampling._select_n_points — addresses the
    `uniformcausticsampling.py _select_n_points()` checklist item in
    issue #65. The method partitions a requested n_points across 1, 2,
    or 3 caustics in proportion to caustic "area", enforcing a per-
    caustic minimum.
    """
    @classmethod
    def setUpClass(cls):
        # __init__ runs a full numerical integration, so we build each
        # topology once and reuse the instance across tests.
        cls.s1 = mm.UniformCausticSampling(s=1.1, q=0.1, n_points=500)
        cls.s2 = mm.UniformCausticSampling(s=2.0, q=0.01, n_points=500)
        cls.s3 = mm.UniformCausticSampling(s=0.5, q=0.001, n_points=500)

    def test_one_caustic_returns_full_n_points(self):
        self.assertEqual(self.s1._n_caustics, 1)
        self.assertEqual(self.s1._select_n_points(100, 10), [100])

    def test_two_caustics_proportional_split_with_min_bumping_n1(self):
        """
        _which_caustic ~ [0, 0.139, 1] -> area_1 << area_2.
        Raw n_1 = 3 falls below min=5 -> bumped to 5, n_2 = 95.
        """
        self.assertEqual(self.s2._n_caustics, 2)
        self.assertEqual(self.s2._select_n_points(100, 5), [5, 95])

    def test_two_caustics_n2_below_min_is_bumped(self):
        """
        Force the n_2 < min branch by patching _which_caustic so area_1
        dominates. Raw n_2 = 0 -> bumped to 10, n_1 = 90.
        """
        original = self.s2._which_caustic
        try:
            self.s2._which_caustic = [0., 0.95, 1.0]
            self.assertEqual(self.s2._select_n_points(100, 10), [90, 10])
        finally:
            self.s2._which_caustic = original

    def test_three_caustics_proportional_split(self):
        """
        _which_caustic ~ [0, 0.132, 0.566, 1].
        fraction_2 ~ 0.477 -> n_2 = 143, n_3 = 143, n_1 = 14.
        """
        self.assertEqual(self.s3._n_caustics, 3)
        self.assertEqual(self.s3._select_n_points(300, 5), [14, 143, 143])

    def test_three_caustics_n1_min_constraint_branch(self):
        """
        When the raw n_1 falls below n_min, n_1 is bumped to n_min and
        n_2/n_3 are recomputed from the remainder as (n-n_1)//2 and
        (n-n_1-n_2). Note: n_2 and n_3 are NOT re-checked against
        n_min, so this branch can produce a result where they fall
        below n_min (here: [11, 9, 10] with n_min=11). Tested as-is
        for coverage; the asymmetry is flagged as a potential
        follow-up in the PR body.
        """
        result = self.s3._select_n_points(30, 11)
        self.assertEqual(result, [11, 9, 10])
        self.assertEqual(len(result), 3)
        self.assertEqual(sum(result), 30)
