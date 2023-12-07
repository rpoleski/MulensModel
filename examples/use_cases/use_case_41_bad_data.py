"""
Use case for changing the bad data array.
"""
import matplotlib.pyplot as plt
import numpy as np

import MulensModel as mm


raise NotImplementedError(
    "We don't know how to enable this functionality.")

data = mm.MulensData(file_name="my_data.dat")
data.plot(show_bad=True)

data.bad[12] = True
data.plot(show_bad=True)

data.bad[np.where(data.mag < 13.)] = True
data.plot(show_bad=True)

data.good[data.mag < 12.] = True
data.plot(show_bad=True)

plt.show()
