from MulensModel.lens import Lens
import unittest
from astropy import units as u

### Test mass setters and getters
def test_2_masses_success():
    lens = Lens()
    lens.mass_1 = 1.0*u.solMass
    lens.mass_2 = 0.3*u.solMass
    assert lens.mass_1 == 1.0*u.solMass
    assert lens.mass_2 == 0.3*u.solMass
    assert lens.total_mass == 1.3*u.solMass

def test_single_mass_success():
    lens_single = Lens()
    lens_single.mass = 0.5*u.solMass
    assert lens_single.mass == 0.5*u.solMass
    assert lens_single.mass_1 == 0.5*u.solMass
    assert lens_single.total_mass == 0.5*u.solMass

def test_distance_success():
    lens = Lens()
    lens.distance = 5.*1000.*u.pc
    assert lens.distance == 5000.*u.pc

class TestQRoutines(unittest.TestCase):
    def test_q_success(self):
        lens = Lens()
        lens.mass_1 = 1.0*u.solMass
        lens.q = 10.**3
        assert lens.mass_2 == 10.**3*u.solMass
        assert lens.total_mass == (1.0+10.**3)*u.solMass
        
class TestInitRoutines(unittest.TestCase):
    def test_init_success(self):
        lens_1 = Lens(mass_1=1.8*u.solMass, mass_2=3.*u.solMass)
        self.assertAlmostEqual(lens_1.mass_1.value, 1.8, places=7)
        self.assertAlmostEqual(lens_1.mass_2.value, 3., places=7)
        
        lens_2 = Lens(mass=1.8*u.solMass, distance=3.*10.**3*u.pc)
        self.assertAlmostEqual(lens_2.mass.value, 1.8, places=7)
        self.assertAlmostEqual(lens_2.mass_1.value, 1.8, places=7)
        self.assertAlmostEqual(lens_2.distance.value, 3.*10.**3,places=7)
        
        lens_3 = Lens(mass_1=1.0*u.solMass, q=0.1,s=0.9)
        self.assertAlmostEqual(lens_3.mass_1.value, 1.0, places=7)
        self.assertAlmostEqual(lens_3.mass_2.value, 0.1, places=7)
        assert lens_3.s == 0.9

        lens_4 = Lens(s=1.0, q=0.1)
        assert lens_4.q == 0.1
        assert lens_4.s == 1.0

def test_a_proj_success():
    lens = Lens(mass_1=1.0*u.solMass,mass_2=0.1*u.solMass,a_proj=1.0*u.au,
                distance=6.*u.kpc)
    assert lens.total_mass == 1.1*u.solMass
    assert lens.q == 0.1

class TestEpsilonRoutines(unittest.TestCase):
    def test_3_body_success(self):
        lens_1 = Lens(total_mass=1.3*u.solMass,q=[0.1,0.2],s=[0.1,0.3])
        
        self.assertAlmostEqual(lens_1.mass_1.value, 1., places=7)
        self.assertAlmostEqual(lens_1.mass_2.value, 0.1, places=7)
        self.assertAlmostEqual(lens_1.mass_3.value, 0.2, places=7)

        lens_2 = Lens()
        lens_2.mass_1 = 0.7*u.solMass
        lens_2.mass_2 = 0.1*u.solMass
        lens_2.mass_3 = 0.2*u.solMass
        
        self.assertAlmostEqual(lens_2.mass_1.value, 0.7, places=7)
        self.assertAlmostEqual(lens_2.mass_2.value, 0.1, places=7)
        self.assertAlmostEqual(lens_2.mass_3.value, 0.2, places=7)
        self.assertAlmostEqual(lens_2.total_mass.value, 1., places=7)

