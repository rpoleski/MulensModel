"""
Script for plotting output from benchmarking script.
"""
from matplotlib import pyplot as plt
import matplotlib
import sys
import numpy as np


def str_to_microsec(str_1, str_2):
    """
    Change the 2 provided numbers into a number of microseconds.
    Examples inputs:
        "100", "us"
        "2.76", "ms"
        "3", "s"
    """
    if str_2 == "us":
        k = 1.
    elif str_2 == "ms":
        k = 1.e3
    elif str_2 == "s":
        k = 1.e6
    else:
        raise ValueError('Unrecognized time format: {:}.'.format(str_2))

    value = float(str_1) * k
    return value


def get_ratio(value_1, sigma_1, value_2, sigma_2):
    """calculates value_1/value_2 and its sigma"""
    value = value_1 / value_2
    sigma = value * np.sqrt((sigma_1/value_1)**2+(sigma_2/value_2)**2)
    return (value, sigma)


in_file = sys.argv[1]

n_all = []
text = []
numpy_mean = dict()
numpy_std_dev = dict()
ids = {
    "MM_static": 0,
    "MM_static-NP": 1,
    "pyLIMA_static": 2,
    "MM_piE": 3,
    "MM_piE-NP": 4,
    "pyLIMA_piE": 5}
common_kwargs = {'lw': 3}
plot_kwargs = {
    0: {"fmt": 'rs', "ms": 5, "label": "MM rectilinear", **common_kwargs},
    1: {"fmt": 's', "color": "orange", "ms": 5,
        "label": "MM rectilinear numpy.sum", **common_kwargs},
    2: {"fmt": 'bs', "ms": 5, "label": "pyLIMA rectilinear", **common_kwargs},
    3: {"fmt": 'ro', "ms": 10, "label": "MM parallax", **common_kwargs},
    4: {"fmt": 'o', "color": "orange", "ms": 10,
        "label": "MM parallax numpy.sum", **common_kwargs},
    5: {"fmt": 'bo', "ms": 10, "label": "pyLIMA parallax", **common_kwargs}
    }

results = None

with open(in_file) as in_data:
    for line_ in in_data.readlines():
        line = line_.split()
        if len(line) < 2 or line[1] != "Mean":
            continue
        full_name = line[0][:-1].split("_")
        name = "_".join(full_name[:-1])
        n_name = full_name[-1]
        mean = str_to_microsec(line[5], line[6])
        std_dev = str_to_microsec(line[8], line[9])
        if name == "numpy":
            n_all.append(n_name)
            numpy_mean[n_name] = mean
            numpy_std_dev[n_name] = std_dev
        else:
            if results is None:
                results = np.zeros((len(ids), len(n_all)))
                results_sigma = np.zeros((len(ids), len(n_all)))
            out = get_ratio(mean, std_dev,
                            numpy_mean[n_name], numpy_std_dev[n_name])
            results[ids[name], n_all.index(n_name)] = out[0]
            results_sigma[ids[name], n_all.index(n_name)] = out[1]

x = np.array([float(n) for n in n_all]) * np.array([0.844, 0.885, 0.929])
for i in range(len(ids)):
    plt.errorbar(x, results[i], results_sigma[i],
                 **{**plot_kwargs[i], 'lw': 3})
    x *= np.array([1.07, 1.05, 1.03])

plt.xlim(50., 30000.)
plt.ylim(0.36, 3.46)
plt.xscale('log')
plt.legend(loc='best', fontsize=15)
matplotlib.rcParams.update({'font.size': 15})
plt.xlabel(r"$N_{points}$")
plt.ylabel(r"$T/T_{numpy,rectilinear}$")
plt.show()

