"""
Define the model by providing mass and distance of lens, source
distance, and relative proper motion. Plot the model light curve after
specifying u_0.

"""
import astropy.units as u
import matplotlib.pyplot as plt

from MulensModel.mulensobjects.lens import Lens
from MulensModel.mulensobjects.source import Source
from MulensModel.mulensobjects.mulenssystem import MulensSystem


# Define a Lens star
my_lens = Lens(mass=0.5*u.solMass, distance=6.e3*u.pc)

# Define a Source Star
my_source = Source(distance=8.e3*u.pc)

# Combine them into a lens system
point_lens = MulensSystem(lens=my_lens, source=my_source)
print('The Lens-Source system WITHOUT proper motion specified:')
print(point_lens)

plt.figure()
point_lens.plot_magnification(u_0=0.1)
plt.title('Magnification Curve')

# Give the objects a relative proper motion
point_lens.mu_rel = 4. * u.mas / u.yr
print('------\nThe Lens-Source system WITH proper motion specified:')
print(point_lens)

plt.figure()
point_lens.plot_magnification(u_0=0.1)
plt.title('Magnification Curve w/Proper Motion Defined')

plt.show()
