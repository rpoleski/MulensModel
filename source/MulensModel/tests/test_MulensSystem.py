import numpy as np
from astropy import units as u

from MulensModel.mulenssystem import MulensSystem
from MulensModel.lens import Lens
from MulensModel.source import Source

kappa = 8.14 * u.mas / u.solMass
lens = {'mass': 0.3 * u.solMass, 'dist': 4 * 10**3 * u.pc}
source = {'dist' : 8 * 10**3 * u.pc}
mu_rel = 3. * u.mas / u.yr
pi_rel = (lens['dist'].to(u.mas, equivalencies=u.parallax()) 
          - source['dist'].to(u.mas, equivalencies=u.parallax()))
thetaE = np.sqrt( kappa * lens['mass'] * pi_rel)
tE = thetaE / mu_rel
"""
print(pi_rel)
print(thetaE)
print(tE)
print(
    'pi_rel: {0}, thetaE: {1}, tE: {2}'.format(
        pi_rel, thetaE, tE.to(u.day)))
"""

def test_mulenssystem():
    test_system = MulensSystem(
        lens=Lens(mass=lens['mass'], distance=lens['dist']), 
        source=Source(distance=source['dist']),
        mu_rel=mu_rel)
    
    assert test_system.pi_rel == pi_rel
    assert test_system.theta_E == thetaE
    assert test_system.t_E == tE.to(u.day)


