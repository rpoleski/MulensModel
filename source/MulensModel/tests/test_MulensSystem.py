import numpy as np
from astropy import units as u

from MulensModel.mulensobjects.mulenssystem import MulensSystem
from MulensModel.mulensobjects.lens import Lens
from MulensModel.mulensobjects.source import Source


def test_mulenssystem():
    kappa = 8.144183118384794 * u.mas / u.solMass
    lens = {'mass': 0.3 * u.solMass, 'dist': 4 * 10**3 * u.pc}
    source = {'dist' : 8 * 10**3 * u.pc}
    mu_rel = 3. * u.mas / u.yr
    pi_rel = (lens['dist'].to(u.mas, equivalencies=u.parallax()) 
              - source['dist'].to(u.mas, equivalencies=u.parallax()))
    thetaE = np.sqrt( kappa * lens['mass'] * pi_rel)
    tE = thetaE / mu_rel

    test_system = MulensSystem(
        lens=Lens(mass=lens['mass'], distance=lens['dist']), 
        source=Source(distance=source['dist']),
        mu_rel=mu_rel)
    
    assert test_system.pi_rel == pi_rel
    assert test_system.theta_E == thetaE
    assert test_system.t_E == tE.to(u.day)

def test_mulenssytem():
    #This test fails due to a discrepancy in the value of kappa (see
    #above value) and kappa as calculated from astropy constants in
    #MulensSystem. Note: 8/23 changed the value of kappa used here and
    #updated theta_E but not subsequent values (test still fails for
    #thetaE).

    lens = Lens(mass=0.64, distance=4.0)
    source = Source(distance=8.0)
    system = MulensSystem(lens=lens, source=source)

    np.testing.assert_almost_equal(system.theta_E.value, 0.807177, decimal=5)
    np.testing.assert_almost_equal(system.r_E.value, 3.22909, decimal=5)
    np.testing.assert_almost_equal(system.r_E_tilde.value, 6.45818, decimal=5)
