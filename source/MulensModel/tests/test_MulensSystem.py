import numpy as np
from astropy import units as u

import MulensModel as mm


def test_mulenssystem_1():
    """test some basic calculations"""
    kappa = 8.144183118384794 * u.mas / u.solMass
    lens = {'mass': 0.3 * u.solMass, 'dist': 4 * 10**3 * u.pc}
    source = {'dist': 8 * 10**3 * u.pc}
    mu_rel = 3. * u.mas / u.yr
    pi_rel = (lens['dist'].to(u.mas, equivalencies=u.parallax())
              - source['dist'].to(u.mas, equivalencies=u.parallax()))
    thetaE = np.sqrt(kappa * lens['mass'] * pi_rel)
    tE = thetaE / mu_rel
    pi_E = pi_rel / thetaE

    test_system = mm.MulensSystem(
        lens=mm.Lens(mass=lens['mass'], distance=lens['dist']),
        source=mm.Source(distance=source['dist']),
        mu_rel=mu_rel)

    assert test_system.pi_rel == pi_rel
    assert abs(test_system.theta_E / thetaE - 1.) < 1.2e-4
    assert abs(test_system.pi_E / pi_E - 1.) < 1.2e-4
    assert isinstance(test_system.pi_E, float)
    assert abs(test_system.t_E / tE - 1.) < 1.2e-4


def test_mulenssystem_2():
    """test basic calculations of theta_E etc."""
    lens = mm.Lens(mass=0.64, distance=4.0)
    source = mm.Source(distance=8.0)
    system = mm.MulensSystem(lens=lens, source=source)

    assert abs(system.theta_E.value / 0.807177 - 1.) < 1.2e-4
    assert abs(system.r_E.value / 3.228708 - 1.) < 1.2e-4
    assert abs(system.r_E_tilde.value / 6.457416 - 1.) < 1.2e-4
