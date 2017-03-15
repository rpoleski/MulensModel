import numpy as np
from astropy import units as u

from MulensModel.mulensobjects.mulenssystem import MulensSystem
from MulensModel.mulensobjects.lens import Lens
from MulensModel.mulensobjects.source import Source

kappa = 8.14611833456212 * u.mas / u.solMass
lens = {'mass': 0.3 * u.solMass, 'dist': 4 * 10**3 * u.pc}
source = {'dist' : 8 * 10**3 * u.pc}
mu_rel = 3. * u.mas / u.yr
pi_rel = (lens['dist'].to(u.mas, equivalencies=u.parallax()) 
          - source['dist'].to(u.mas, equivalencies=u.parallax()))
thetaE = np.sqrt( kappa * lens['mass'] * pi_rel)
tE = thetaE / mu_rel

def test_mulenssystem():
    test_system = MulensSystem(
        lens=Lens(mass=lens['mass'], distance=lens['dist']), 
        source=Source(distance=source['dist']),
        mu_rel=mu_rel)
    
    assert test_system.pi_rel == pi_rel
    assert test_system.theta_E == thetaE
    assert test_system.t_E == tE.to(u.day)

def test_mulenssytem():
    lens = Lens(mass=0.64, distance=4.0)
    source = Source(distance=8.0)
    system = MulensSystem(lens=lens, source=source)

    np.testing.assert_almost_equal(system.theta_E.value, 0.80727, decimal=5)
    np.testing.assert_almost_equal(system.r_E.value, 3.22909, decimal=5)
    np.testing.assert_almost_equal(system.r_E_tilde.value, 6.45818, decimal=5)
