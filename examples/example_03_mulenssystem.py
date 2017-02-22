import astropy.units as u
import matplotlib.pyplot as pl

from MulensModel.lens import Lens
from MulensModel.source import Source
from MulensModel.mulenssystem import MulensSystem

# Define a Lens star
my_lens = Lens(mass=0.5*u.solMass, distance=6.e3*u.pc)

# Define a Source Star
my_source = Source(distance=8.e3*u.pc)

# Combine them into a lens system
point_lens = MulensSystem(lens=my_lens, source=my_source)
print(point_lens)
pl.figure()
point_lens.plot_magnification(u_0=0.1)
pl.title('Magnification Curve')

# Given the objects a relative proper motion
point_lens.mu_rel = 4. * u.mas / u.yr
print(point_lens)
pl.figure()
point_lens.plot_magnification(u_0=0.1)
pl.title('Magnification Curve w/Proper Motion Defined')

pl.show()
