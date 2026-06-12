import os
import unittest

import numpy as np

import MulensModel as mm


FILE_01 = os.path.join(mm.DATA_PATH, "unit_test_files",
                       "MB11293_caustics.dat")
test_caustics = np.genfromtxt(FILE_01, names=['X', 'Y'], dtype=None)


def get_index_of_nearest_point(test_data, x, y):
    """
    Find index of test_data (*dict*) that is closes to (x, y)
    """
    return np.argmin((test_data['X']-x)**2+(test_data['Y']-y)**2)


def test_caustic():
    """
    Make sure the caustic is properly calculated.
    """
    s = 0.548
    q = 0.0053

    caustics = mm.CausticsBinary(q=q, s=s)

    x, y = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = get_index_of_nearest_point(test_caustics, x[i], y[i])

        np.testing.assert_almost_equal(
            x[i], test_caustics['X'][index], decimal=5)
        np.testing.assert_almost_equal(
            y[i], test_caustics['Y'][index], decimal=5)


class TestCriticalCurve(unittest.TestCase):
    """
    Tests for CausticsBinary.critical_curve — addresses the
    `caustics.py critical_curve()` checklist item in issue #65.
    """
    def setUp(self):
        self.caustics = mm.CausticsBinary(q=0.0053, s=0.548)

    def test_critical_curve_returns_object_with_x_y_lists(self):
        cc = self.caustics.critical_curve
        self.assertIsInstance(cc.x, list)
        self.assertIsInstance(cc.y, list)
        self.assertEqual(len(cc.x), len(cc.y))
        # n_points=5000 default -> n_angles ~1250, 4 roots each
        self.assertEqual(len(cc.x), 5000)

    def test_critical_curve_consistent_across_calls(self):
        """Two accesses return the same x, y values."""
        first = self.caustics.critical_curve
        second = self.caustics.critical_curve
        self.assertEqual(first.x, second.x)
        self.assertEqual(first.y, second.y)


class TestCausticsBinaryInit(unittest.TestCase):
    """
    Constructor edge cases for CausticsBinary.q (same area as #65 item
    above; covers the q-as-list branch in __init__).
    """
    def test_q_as_single_element_list_is_accepted(self):
        caustics = mm.CausticsBinary(q=[0.0053], s=0.548)
        self.assertEqual(caustics.q, 0.0053)

    def test_q_as_multi_element_list_raises_not_implemented(self):
        with self.assertRaisesRegex(NotImplementedError, "Too many q"):
            mm.CausticsBinary(q=[0.0053, 0.1], s=0.548)
