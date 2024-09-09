import numpy as np
import unittest

import MulensModel as mm


#def test_small_q():
#    """
#    run calculations with small mass-ratio
#    """
#    q = 1.e-8
#    s = 1.8
#    x = s - 1. / s
#    y = 5.e-7 * 0.
#
#    trajectory = make_trajectory((x, y), {'s': s, 'q': q})
#    lens = mm.BinaryLensPointSourceWM95Magnification(trajectory=trajectory)
#    result = lens.get_magnification()
#    np.testing.assert_almost_equal(result, 3.6868957, decimal=3)
#
#
#def test_binary_lens_hexadecapole():
#    """
#    tests hexadecapole and quadrupole calculation for planetary case
#    with rho=0.001 and 3 values of gamma: 0, 0.5, and 1.0
#    """
#    s = 1.35
#    q = 0.00578
#    rho = 0.001
#    x = 0.7142010570568691 - s * q / (1. + q)
#    y = 0.00189679191923936
#
#    trajectory = make_trajectory((x, y), {'s': s, 'q': q, 'rho': rho})
#
#    pspl_mag = 4.691830779895085
#    reference_00 = [5.017252440557196, 4.9587638949353146, pspl_mag]
#    reference_05 = [4.981368071884021, 4.932070583431292, pspl_mag]
#    reference_10 = [4.9454837032108445, 4.905377271927269, pspl_mag]
#    references = [reference_00, reference_05, reference_10]
#
#    gammas = [0.0, 0.5, 1.0]
#    for i in range(3):
#        print(i)
#        bl = mm.BinaryLensHexadecapoleMagnification(
#            trajectory=trajectory, gamma=gammas[i],
#            all_approximations=True)
#
#        result = bl.get_magnification()
#        np.testing.assert_almost_equal(result, references[i])
#
#
#def test_vbbl_1():
#    """
#    check basic magnification calculation using VBBL
#    """
#    s = 0.8
#    q = 0.1
#
#    m_1 = 1. / (1. + q)
#    m_2 = q / (1. + q)
#    bl = mm.BinaryLens(m_1, m_2, s)
#
#    result = bl.vbbl_magnification(0.01, 0.01, 0.01)
#    np.testing.assert_almost_equal(result, 18.2834436, decimal=3)
#
#
#def test_vbbl_2():
#    """
#    Check VBBL magnification calculation for binary lens with finite source
#    that was producing wrong results for VBBL earlier than v3.5
#    """
#    x = [-2.8798499936424813, -2.87980198609534, -2.879750341503788]
#    y = [0.2603315602357186, 0.26034667859291694, 0.26036294250727565]
#    s = 0.3121409537799967
#    q = 0.0018654668855723224
#    rho = 0.002966662955047919
#
#    m_1 = 1. / (1. + q)
#    m_2 = q / (1. + q)
#    bl = mm.BinaryLens(m_1, m_2, s)
#    results = [bl.vbbl_magnification(x_, y_, rho) for (x_, y_) in zip(x, y)]
#    # VBBL 2.0.1 was returning:
#    # [1.345365452870409, 1.368843518228974, 1.3442156685350604]
#    # i.e., the second value was wrong.
#    # VBBL 3.5 returns:
#    # [1.3455151798453464, 1.3449680490180915, 1.344226470628085]
#    np.testing.assert_almost_equal(results[2], 1.344226470628085, decimal=5)
#    mean = (results[0] + results[2]) / 2
#    np.testing.assert_almost_equal(results[1], mean, decimal=4)
#
#
#def test_ac():
#    """
#    check basic magnification calculation using AdaptiveContouring
#    """
#    s = 0.8
#    q = 0.1
#
#    m_1 = 1. / (1. + q)
#    m_2 = q / (1. + q)
#    bl = mm.BinaryLens(m_1, m_2, s)
#
#    # Previous version was failing because of bug in AC:
#    # result = bl.adaptive_contouring_magnification(
#    #     0.01, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
#    # np.testing.assert_almost_equal(result, 18.2834436, decimal=3)
#    result = bl.adaptive_contouring_magnification(
#        0.06, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
#    np.testing.assert_almost_equal(result, 11.403036510555962, decimal=3)


# Magnification Object Tests

def make_trajectory(pos, params):
    times = pos[0]
    u0 = pos[1]

    t0 = 0
    tE = 1
    alpha = 180.
    ulens_params = {
        't_0': t0, 'u_0': u0, 't_E': tE,
        's': params['s'], 'q': params['q'], 'alpha': alpha}
    if 'rho' in params:
        ulens_params['rho'] = params['rho']

    trajectory = mm.Trajectory(
        times=times, parameters=mm.ModelParameters(ulens_params))

    return trajectory


class TestBinaryLensPointSourceMagnification(unittest.TestCase):

    def setUp(self):
        """
        run calculations with small mass-ratio
        """
        self.q = 1.e-8
        self.s = 1.8
        self.x = self.s - 1. / self.s
        self.y = 5.e-7 * 0.
        self.expected = 3.6868957

        self.trajectory = make_trajectory(
            [self.x, self.y], {'s': self.s, 'q': self.q})

    def test_BLPS_WittMao95(self):
        lens = mm.BinaryLensPointSourceWM95Magnification(
            trajectory=self.trajectory)
        result = lens.get_magnification()
        np.testing.assert_almost_equal(result, self.expected, decimal=3)

    def test_BLPS_VBBL(self):
        lens = mm.BinaryLensPointSourceVBBLMagnification(
            trajectory=self.trajectory)
        result = lens.get_magnification()
        np.testing.assert_almost_equal(result, self.expected, decimal=3)


def test_BinaryLensQuadrupoleMagnification():
    raise NotImplementedError('No tests defined for quadrupole approximation!')
    assert 1 == 2


class TestBinaryLensHexadecapoleMagnification(unittest.TestCase):
    """
    tests hexadecapole and quadrupole calculation for planetary case
    with rho=0.001 and 3 values of gamma: 0, 0.5, and 1.0
    """

    def setUp(self):
        self.s = 1.35
        self.q = 0.00578
        self.rho = 0.001
        self.x = 0.7142010570568691 - self.s * self.q / (1. + self.q)
        self.y = 0.00189679191923936

        self.trajectory = make_trajectory(
            [self.x, self.y], {'s': self.s, 'q': self.q, 'rho': self.rho})

        self.pspl_mag = 4.691830779895085
        # [hexa, quad, pspl]
        self.reference_00 = [5.017252440557196, 4.9587638949353146, self.pspl_mag]
        self.reference_05 = [4.981368071884021, 4.932070583431292, self.pspl_mag]
        self.reference_10 = [4.9454837032108445, 4.905377271927269, self.pspl_mag]

    def _test_gamma(self, gamma, reference):
        """
        Check if hexadecapole calculation with all_approximations works well.
        """
        lens = mm.BinaryLensHexadecapoleMagnification(
            trajectory=self.trajectory, gamma=gamma, all_approximations=True)
        result = np.array(lens.get_magnification()).flatten()
        np.testing.assert_almost_equal(result, reference)

    def test_gamma_00(self):
        self._test_gamma(0.0, self.reference_00)

    def test_gamma_05(self):
        self._test_gamma(0.5, self.reference_05)

    def test_gamma_10(self):
        self._test_gamma(1.0, self.reference_10)


class TestBinaryLensQuadrupoleMagnification(
    TestBinaryLensHexadecapoleMagnification):

    def setUp(self):
        TestBinaryLensHexadecapoleMagnification.setUp(self)

    def test_gamma_00(self):
        print('gamma00', self.trajectory)
        lens = mm.BinaryLensQuadrupoleMagnification(
            trajectory=self.trajectory, gamma=0.0)
        result = lens.get_magnification()
        np.testing.assert_almost_equal(result, self.reference_00[1])

    def test_gamma_05(self):
        pass

    def test_gamma_10(self):
        pass

def test_BinaryLensVBBLMagnification_1():
    """
    check basic magnification calculation using VBBL
    """
    s = 0.8
    q = 0.1
    x = 0.01
    y = 0.01
    rho = 0.01

    trajectory = make_trajectory([x, y], {'s': s, 'q': q, 'rho': rho})
    lens = mm.BinaryLensVBBLMagnification(trajectory=trajectory)
    result = lens.get_magnification()
    np.testing.assert_almost_equal(result, 18.2834436, decimal=3)


def test_BinaryLensVBBLMagnification_2():
    """
    Check VBBL magnification calculation for binary lens with finite source
    that was producing wrong results for VBBL earlier than v3.5
    """
    x = [-2.8798499936424813, -2.87980198609534, -2.879750341503788]
    y = [0.2603315602357186, 0.26034667859291694, 0.26036294250727565]
    s = 0.3121409537799967
    q = 0.0018654668855723224
    rho = 0.002966662955047919

    results = []
    for (x_, y_) in zip(x, y):
        trajectory = make_trajectory([x_, y_], {'s': s, 'q': q, 'rho': rho})
        lens = mm.BinaryLensVBBLMagnification(trajectory=trajectory)
        result = lens.get_magnification()
        results.append(result)

    # VBBL 2.0.1 was returning:
    # [1.345365452870409, 1.368843518228974, 1.3442156685350604]
    # i.e., the second value was wrong.
    # VBBL 3.5 returns:
    # [1.3455151798453464, 1.3449680490180915, 1.344226470628085]
    np.testing.assert_almost_equal(results[2], 1.344226470628085, decimal=5)
    mean = (results[0] + results[2]) / 2
    np.testing.assert_almost_equal(results[1], mean, decimal=4)


def test_BinaryLensVBBLMagnification_LD():
    raise NotImplementedError('No tests defined for VBBL + gamma!')
    assert 1 == 2


def test_BinaryLensAdaptiveContouringMagnification():
    """
    check basic magnification calculation using AdaptiveContouring
    """
    s = 0.8
    q = 0.1
    x = 0.06
    y = 0.01
    rho = 0.01

    # Previous version was failing because of bug in AC:
    # result = bl.adaptive_contouring_magnification(
    #     0.01, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
    # np.testing.assert_almost_equal(result, 18.2834436, decimal=3)
    trajectory = make_trajectory([x, y], {'s': s, 'q': q, 'rho': rho})
    lens = mm.BinaryLensAdaptiveContouringMagnification(
        trajectory=trajectory,  accuracy=0.019, ld_accuracy=1e-3)
    result = lens.get_magnification()
    np.testing.assert_almost_equal(result, 11.403036510555962, decimal=3)


if __name__ == '__main__':
    test = TestBinaryLensPointSourceMagnification()
    test.setUp()
    test.test_BLPS_WittMao95()
