"""
Show plotting caustics using different methods.
"""
import matplotlib.pyplot as plt
import numpy as np

import MulensModel as mm


s = 1.2
q = 0.5
n_points = 200

color = np.linspace(0., 1., n_points)

# First, we use standard procedure to plot the caustic. You will see that
# the points are distributed non-uniformly, i.e., the density is higher
# near cusps. We also add color to indicate order of plotting.
# It turns out to be a complicated shape.
caustic = mm.CausticsBinary(s=s, q=q)
caustic.plot(c=color, n_points=n_points)
plt.axis('equal')
plt.colorbar()
plt.title('standard plotting')

# Second, we use uniform sampling. The density of points is constant.
# Here the color scale indicates x_caustic values.
plt.figure()
sampling = mm.UniformCausticSampling(s=s, q=q)
points = [sampling.caustic_point(c) for c in color]
x = [p.real for p in points]
y = [p.imag for p in points]
plt.scatter(x, y, c=color)
plt.axis('equal')
plt.colorbar()
plt.title('uniform sampling plotting')

plt.show()
