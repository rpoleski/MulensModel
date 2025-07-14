"""
Script to prepare files with prior used as input `from file` in ulens_model_fit.py base on simulated microlensing
events using genulens code (https: //github.com/nkoshimoto/genulens)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rv_histogram   # ,gaussian_kde
import fastkde   # https: // github.com/LBL-EESA/fastkde.git


def Gauss(x, A, m, s):
    return A*(1./s/np.sqrt(2*np.pi)) * np.exp(-((x-m)**2) / (2 * s**2))


def get_data(file):
    """
    Read data from file, skipping header and return it as a structured numpy array.
    """
    header_line = None
    with open(file, "r") as f:
        for line in f:
            if line.strip().startswith("#"):
                header_line = line
            else:
                break
    # Extract column names from the header line
    column_names = header_line.lstrip("#").strip().split()
    print("Column names:", column_names)
    data = np.genfromtxt(file, comments='#', names=column_names, dtype=np.float64)
    return data


def get_pdf(values, weights, parameter, nbins=10000):
    """
    Function to get the PDF from the simulated values.
    """
    histogram, bin_edges = np.histogram(values, bins=nbins, weights=weights)
    rv = rv_histogram((histogram, bin_edges))
    pdf_orginal = rv.pdf
    dbins = (bin_edges[1] - bin_edges[0]) / 2
    bins = (bin_edges+dbins)[:-1]
    histogram = spread(bins, histogram, parameter)
    spreaded = rv_histogram([histogram, bin_edges])
    limits = spreaded.interval(0.98)
    pdf_spread_smooth = smooth(spreaded.rvs(size=10000000), bins)
    return pdf_orginal, pdf_spread_smooth, bins, limits


def spread(values, counts, parameter):
    """
    Function to spread tails of the simulated distributions for better sampling.
    It fits a gaussian function to the histogram, multiplies its sigma by a scale factor,
    and adds the gaussian function to the initial histogram.
    """
    scale = 4.
    mean0 = np.average(values, weights=counts)
    sigma0 = np.sqrt(np.average(
        (values - mean0)**2, weights=counts))
    p0 = [max(counts), mean0, sigma0]
    pars, _ = curve_fit(Gauss, values, counts, p0=p0)
    fit_y = Gauss(values, pars[0], pars[1], pars[2] * scale)
    print('{:s}:  A= {:.2f} mean= {:.2f}, sigma={:.2f}'.format(
        parameter, pars[0], pars[1], pars[2]*scale))
    spread_histogram = np.add(counts, fit_y)
    return spread_histogram


def smooth(sample, bins):
    """
    Function that estimate PDF based on the sample using a gaussian kernel density estimation.
    It returns the smoothed PDF.
    """
    pdf = fastkde.fastKDE.pdf_at_points(sample, list_of_points=bins)

    def pdf_spread_smooth(x):
        return np.interp(x, bins, pdf)
    # slower way and less accurate:
    # smooth = gaussian_kde(sample)
    # pdf_spread_smooth = smooth.pdf
    return pdf_spread_smooth


def plot_prior(name, limits,  pdf_orginal, pdf_spread_smooth, prefix_output):
    true = np.linspace(limits[0], limits[1], 10000)
    plt.plot(true, pdf_orginal(true), label='PDF friom genulens')
    plt.plot(true, pdf_spread_smooth(true), label='PDF spreaded and smoothed')
    plt.title('PDF '+name)
    plt.legend()
    file = prefix_output+name+'.png'
    print("Saving plot of used prior on " + name + " to " + file)
    plt.savefig(file)
    plt.close()


simulated_file = 'data/OB03235/OB03235_genulens.out'
prefix_output = 'data/OB03235/OB03235_prior_'
data = get_data(simulated_file)
weights = data['wtj']

# choose parameters
parameters = ['t_E', 'pi_E_E', 'pi_E_N']

# convesion from MulensModel to genulens names
conve = {'pi_E_E': 'pi_EE',
         'pi_E_N': 'pi_EN',
         }

for parameter in parameters:
    parameter_genulens = conve.get(parameter, parameter)
    values = data[parameter_genulens]
    pdf_orginal, pdf_spread_smooth, bins, limits = get_pdf(values, weights, parameter)
    plot_prior(parameter, limits, pdf_orginal, pdf_spread_smooth, prefix_output)
    file_out = prefix_output + parameter + '.txt'
    combined = np.column_stack((bins, pdf_spread_smooth(bins)))
    np.savetxt(file_out, combined, fmt='%.6f', delimiter=' ', header='# x pdf(x)')
