import numpy as np

from MulensModel import Model

def test_orbital_motion():
    """test orbital motion parameters"""
    params = {'t_0': 2456789.01234, 'u_0': 0.01, 't_E': 25., 'alpha': 30., 
        's': 1.2345, 'q': 0.555}
    model_static = Model(params)
    model_motion = Model({**params, 'dalpha_dt': 10., 
        'ds_dt': 0.1})

    # Test is_static()
    assert model_static.is_static()
    assert not model_motion.is_static()
    assert model_static.parameters.is_static()
    assert not model_motion.parameters.is_static()

    epoch_1 = params['t_0'] - 18.2625
    epoch_2 = params['t_0']
    epoch_3 = params['t_0'] + 18.2625

    # Test get_s() for static case
    static = model_static.parameters
    assert static.get_s(epoch_1) == static.get_s(epoch_2)
    assert static.get_s(epoch_1) == static.get_s(epoch_3)
    assert static.get_s(epoch_1) == params['s']

    # Test get_alpha() for orbital motion case
    motion = model_motion.parameters
    np.testing.assert_almost_equal(motion.get_alpha(epoch_1).value, 29.5)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_2).value, 30.)
    np.testing.assert_almost_equal(motion.get_alpha(epoch_3).value, 30.5)

    # Test get_s() for orbital motion case
    np.testing.assert_almost_equal(motion.get_s(epoch_1), 1.2295)
    np.testing.assert_almost_equal(motion.get_s(epoch_2), 1.2345)
    np.testing.assert_almost_equal(motion.get_s(epoch_3), 1.2395)

# We further need to test .gamma_parallel .gamma_perp .gamma,
# and causitcs for epoch=XXX,
# and some version of t_0_orb.

def test_t_binary():
    """Test that t_binary is correctly set"""
    # If not defined, t_binary = t_0
    pass

def test_magnification():
    """Test that the magnification is correcly calculated for an orbital
    motion case"""
    pass
