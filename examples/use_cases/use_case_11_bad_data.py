import matplotlib.pyplot as plt

import MulensModel as mm


# Needs commentary

raise NotImplementedError(
    'Most of this use case is not implemented. Need to reconsider ' +
    'estimate_model_params() given that Model cannot be defined this way.')

data = mm.MulensData(file_name="my_data.dat")
model = mm.Model(n_components=1)  # This is not allowed

data.bad[np.isnan(data.err)] = True

event = mm.Event(datasets=data, model=model)
event.estimate_model_params()  # aspirational
event.get_chi2()
event.clean_data()
event.clean_data(sigma=3.)
"""
if sigma==None, set sigma based on the number of data points.
clean_data requires fitting a model...
"""
event.get_chi2()  # should mask bad data

plt.scatter(data.time, data.mag, marker="o", facecolor=None)
plt.scatter(data.time[data.good], data.mag[not data.bad],
            marker="o", facecolor="black")
plt.scatter(data.time[data.bad], data.mag[data.bad], marker="x")
