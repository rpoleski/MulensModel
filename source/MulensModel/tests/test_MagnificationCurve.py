import numpy as np

from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.modelparameters import ModelParameters


def test_fspl_noLD():
    """check if FSPL magnification is calculate properly"""
    t_0 = 2456789.012345
    t_E = 23.4567
    u_0 = 1e-4
    rho = 1e-3
    t_vec = np.array([-(rho**2-u_0**2)**0.5, 0., ((0.5*rho)**2-u_0**2)**0.5])
    t_vec = t_vec * t_E + t_0
    
    params = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_E':t_E, 'rho':rho})

    mag_curve = MagnificationCurve(times=t_vec, parameters=params)
    methods = [t_0-t_E, 'finite_source_uniform_Gould94', t_0+t_E]
    mag_curve.set_magnification_methods(methods, 'point_source')
    results = mag_curve.get_point_lens_magnification()
    
    u = np.array([rho, u_0, 0.5*rho])
    pspl = (u**2 + 2.) / np.sqrt(u**2 * (u**2 + 4.))
    expected = np.array([1.27323965, 0.19949906, 0.93421546]) # These values 
    # were calculated by Andy Gould (file b0b1.dat).
    expected *= pspl
    
    np.testing.assert_almost_equal(expected, results, decimal=4)
    
def test_fspl():
    """check if FSPL magnification is calculate properly"""
    t_0 = 2456789.012345
    t_E = 23.4567
    u_0 = 1e-4
    rho = 1e-3
    gamma = 0.56789
    t_vec = np.array([-(rho**2-u_0**2)**0.5, 0., ((0.5*rho)**2-u_0**2)**0.5])
    t_vec = t_vec * t_E + t_0
    
    params = ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})

    mag_curve = MagnificationCurve(times=t_vec, parameters=params, gamma=gamma)
    methods = [t_0-t_E, 'finite_source_LD_Yoo04', t_0+t_E]
    mag_curve.set_magnification_methods(methods, 'point_source')
    results = mag_curve.get_point_lens_magnification()
    
    u = np.array([rho, u_0, 0.5*rho])
    pspl = (u**2 + 2.) / np.sqrt(u**2 * (u**2 + 4.))
    expected = np.array([1.27323965-gamma*0.09489869, 
                         0.19949906-gamma*-0.03492121, 
                         0.93421546-gamma*-0.09655794]) # These values 
    # were calculated by Andy Gould (file b0b1.dat).
    expected *= pspl
    
    np.testing.assert_almost_equal(expected/results, 1., decimal=4)

def test_Lee09():
    """
    test Lee+2009 finite source calculation
    """
# uniform rho = 1 (u, A):
# 3.53553390593 1.01084060513
# 2.06155281281 1.06962639343
# 1.11803398875 1.42451408166
# 0.707106781187 2.02334097551
# 0.5 2.13919086656

# LD gamma=0.5:
# 3.50142828 1.01110609638
# 2.00249843945 1.07461016241
# 1.00498756211 1.57232954942
# 0.509901951359 2.21990790526
# 0.1 2.39458814753

# 2LD lambda=0.3:
# 3.50142828 1.0110829794
# 2.00249843945 1.07404148634
# 1.00498756211 1.55620547462
# 0.509901951359 2.24809136704
# 0.1 2.44503143812

    params = ModelParameters({'t_0': 0., 'u_0': 0.1, 't_E': 10., 'rho': 0.5})
    t_vec = np.array([0])
    mag_curve = MagnificationCurve(times=t_vec, parameters=params, gamma=0.3)
    methods = [-1., 'finite_source_uniform_Lee09', 1.]
    mag_curve.set_magnification_methods(methods, 'point_source')
    results = mag_curve.get_point_lens_magnification()
    assert results[0] > 4. and results[0] < 4.3

def test_PSPL_for_binary():
    """test PSPL model used in a model that is defined as binary"""
    t_0 = 1000.
    t_E = 20.
    u_0 = 1.
    t_vec = np.array([10., 100.]) * t_E + t_0
    params = ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 's': 1.2, 'q': 0.1, 'alpha': 0.})
    mag_curve = MagnificationCurve(times=t_vec, parameters=params)
    mag_curve.set_magnification_methods(None, 'point_source_point_lens')
    u2 = u_0**2 + ((t_vec - t_0) / t_E)**2
    pspl = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    np.testing.assert_almost_equal(pspl, mag_curve.magnification)

