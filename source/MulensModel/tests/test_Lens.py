import numpy as np
from astropy import units as u

import MulensModel as mm


# Test mass setters and getters
def test_2_masses_success():
    lens = mm.Lens()
    lens.mass_1 = 1.0*u.solMass
    lens.mass_2 = 0.3*u.solMass
    assert lens.mass_1 == 1.0*u.solMass
    assert lens.mass_2 == 0.3*u.solMass
    assert lens.total_mass == 1.3*u.solMass


def test_single_mass_success():
    lens_single = mm.Lens()
    lens_single.mass = 0.5*u.solMass
    assert lens_single.mass == 0.5*u.solMass
    assert lens_single.mass_1 == 0.5*u.solMass
    assert lens_single.total_mass == 0.5*u.solMass


def test_distance_success():
    lens = mm.Lens()
    lens.distance = 5.*1000.*u.pc
    assert lens.distance == 5000.*u.pc


def test_total_mass_q():
    lens = mm.Lens()
    lens.total_mass = 0.8*u.solMass
    lens.q = 0.25
    np.testing.assert_almost_equal(lens.mass_1.value, 0.64)
    np.testing.assert_almost_equal(lens.mass_2.value, 0.16)


def test_q_total_mass():
    lens = mm.Lens()
    lens.q = 0.25
    lens.total_mass = 0.8*u.solMass
    np.testing.assert_almost_equal(lens.mass_1.value, 0.64)
    np.testing.assert_almost_equal(lens.mass_2.value, 0.16)


def test_q_success():
    lens = mm.Lens()
    lens.mass_1 = 1.0*u.solMass
    lens.q = 10.**3
    assert lens.mass_2 == 10.**3*u.solMass
    assert lens.total_mass == (1.0+10.**3)*u.solMass


def test_init_success():
    lens_1 = mm.Lens(mass_1=1.8*u.solMass, mass_2=3.*u.solMass)
    np.testing.assert_almost_equal(lens_1.mass_1.value, 1.8)
    np.testing.assert_almost_equal(lens_1.mass_2.value, 3.)

    lens_2 = mm.Lens(mass=1.8*u.solMass, distance=3.*10.**3*u.pc)
    np.testing.assert_almost_equal(lens_2.mass.value, 1.8)
    np.testing.assert_almost_equal(lens_2.mass_1.value, 1.8)
    np.testing.assert_almost_equal(lens_2.distance.value, 3.*10.**3)

    lens_3 = mm.Lens(mass_1=1.0*u.solMass, q=0.1, s=0.9)
    np.testing.assert_almost_equal(lens_3.mass_1.value, 1.0)
    np.testing.assert_almost_equal(lens_3.mass_2.value, 0.1)
    assert lens_3.s == 0.9

    lens_4 = mm.Lens(s=1.0, q=0.1)
    assert lens_4.q == 0.1
    assert lens_4.s == 1.0


def test_a_proj_success():
    lens = mm.Lens(mass_1=1.0*u.solMass, mass_2=0.1*u.solMass, a_proj=1.0*u.au,
                   distance=6.*u.kpc)
    assert lens.total_mass == 1.1*u.solMass
    assert lens.q == 0.1


def test_3_body_success():
    lens_1 = mm.Lens(
        total_mass=1.3*u.solMass,
        q=[0.1, 0.2],
        s=[0.1, 0.3])

    np.testing.assert_almost_equal(lens_1.mass_1.value, 1.)
    np.testing.assert_almost_equal(lens_1.mass_2.value, 0.1)
    np.testing.assert_almost_equal(lens_1.mass_3.value, 0.2)

    lens_2 = mm.Lens()
    lens_2.mass_1 = 0.7*u.solMass
    lens_2.mass_2 = 0.1*u.solMass
    lens_2.mass_3 = 0.2*u.solMass

    np.testing.assert_almost_equal(lens_2.mass_1.value, 0.7)
    np.testing.assert_almost_equal(lens_2.mass_2.value, 0.1)
    np.testing.assert_almost_equal(lens_2.mass_3.value, 0.2)
    np.testing.assert_almost_equal(lens_2.total_mass.value, 1.)
