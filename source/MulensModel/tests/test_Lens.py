from MulensModel.lens import Lens
import unittest
from astropy import units as u

### Test mass setters and getters
def test_2_masses_success():
    lens = Lens(n_components=2)
    lens.mass_1 = 1.0*u.solMass
    lens.mass_2 = 0.3*u.solMass
    assert lens.mass_1 == 1.0*u.solMass
    assert lens.mass_2 == 0.3*u.solMass
    assert lens.mass == 1.3*u.solMass
    assert lens.total_mass == 1.3*u.solMass

def test_single_mass_success():
    lens_single = Lens(n_components=1)
    lens_single.mass = 0.5*u.solMass
    assert lens_single.mass == 0.5*u.solMass
    assert lens_single.mass_1 == 0.5*u.solMass
    assert lens_single.total_mass == 0.5*u.solMass

def test_distance_success():
    lens = Lens(n_components=1)
    lens.distance = 5.*1000.*u.pc
    assert lens.distance == 5000.*u.pc

class TestQRoutines(unittest.TestCase):
    def test_q_success(self):
        lens = Lens(n_components=2)
        lens.mass_1 = 1.0*u.solMass
        lens.q = 10.**3
        assert lens.mass_2 == 10.**3*u.solMass
        assert lens.mass == (1.0+10.**3)*u.solMass
        
        lens_2 = Lens(n_components=2)
        lens_2.total_mass = 1.0*u.solMass
        lens_2.q = 0.1/0.9
        self.assertAlmostEqual(lens_2.mass_1.value, 0.9, places=7)
        self.assertAlmostEqual(lens_2.mass_2.value, 0.1, places=7)

def test_init_success():
    lens_1 = Lens(mass_1=1.8*u.solMass, mass_2=3.*u.solMass)
    assert lens_1.mass_1 == 1.8*u.solMass
    assert lens_1.mass_2 == 3.*u.solMass
    assert lens_1.n_components == 2

    lens_2 = Lens(mass=1.8*u.solMass, distance=3.*10.**3*u.pc)
    assert lens_2.mass == 1.8*u.solMass
    assert lens_2.mass_1 == 1.8*u.solMass
    assert lens_2.distance == 3.*10.**3*u.pc

    lens_3 = Lens(mass_1=1.0*u.solMass, q=0.1,s=0.9)
    assert lens_3.mass_1 == 1.0*u.solMass
    assert lens_3.mass_2 == 0.1*u.solMass
    assert lens_3.n_components == 2
    assert lens_3.s == 0.9

    lens_4 = Lens(s=1.0, q=0.1)
    assert lens_4.q == 0.1
    assert lens_4.s == 1.0
    assert lens_4.n_components == 2

def test_a_proj_success():
    lens = Lens(mass_1=1.0*u.solMass,mass_2=0.1*u.solMass,a_proj=1.0*u.au,
                distance=6.*u.kpc)
    assert lens.q == 0.1
