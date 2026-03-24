"""
Script to prepare files with prior used as input `from file` in ulens_model_fit.py base on simulated microlensing
events using genulens code (https: //github.com/nkoshimoto/genulens)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import rv_histogram   # ,gaussian_kde
import fastkde   # https: // github.com/LBL-EESA/fastkde.git
from matplotlib.gridspec import GridSpec


def Gauss(x, A, m, s):
    return A*(1./s/np.sqrt(2*np.pi)) * np.exp(-((x-m)**2) / (2 * s**2))


def Gauss_2d(coords, A, mu_x, mu_y, sigma_x, sigma_y):
    x, y = coords
    return A * np.exp(-(((x - mu_x)**2) / (2 * sigma_x**2) + ((y - mu_y)**2) / (2 * sigma_y**2)))


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


def get_pdf_2D(values_x, values_y, weights, parameter_2D, prefix_output, nbins=1000):
    """
    Function to get the PDF from the simulated values.
    """
    histogram, xedges, yedges = np.histogram2d(values_x, values_y, bins=nbins, weights=weights)
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(xcenters, ycenters)
    histogram = spread_2D(X, Y, histogram, parameter_2D)
    values_spread_x, values_spread_y = sample_2D(histogram, nbins, xedges, yedges)
    plot_2D_points(values_x, values_y, values_spread_x, values_spread_y, parameter_2D, prefix_output)
    pdf_spread_smooth, x, y = smooth_2D(values_spread_x, values_spread_y, nbins=nbins)
    safe_pdf_2D(pdf_spread_smooth, x, y, parameter_2D, prefix_output)
    plot_pdf_2D(pdf_spread_smooth, x, y, parameter_2D, prefix_output)
    return pdf_spread_smooth, xcenters, ycenters, histogram


def plot_2D_points(values_x, values_y, values_spread_x, values_spread_y, parameter_2D, prefix_output):
    plt.scatter(values_spread_x, values_spread_y, alpha=0.5, s=0.01, label='Spreaded')
    plt.scatter(values_x, values_y, alpha=0.5, s=0.01, label='Original')

    plt.xlabel(parameter_2D[0])
    plt.ylabel(parameter_2D[1])
    plt.title('Scatter plot of original and spreaded samples')
    plt.legend()
    file = prefix_output + 'scatter_' + parameter_2D[0] + '_' + parameter_2D[1] + '.png'
    print("Saving scatter plot of original and spreaded samples to " + file)
    plt.savefig(file)
    plt.close()


def safe_pdf_2D(pdf_spread_smooth, x, y, parameter_2D, prefix_output):
    z = pdf_spread_smooth
    if len(z) != len(x) and len(z) != len(y):
        raise ValueError("Length of z does not match length of x and y")
    data = np.column_stack((x, y, z))
    file_out = prefix_output + parameter_2D[0] + '_' + parameter_2D[1] + '.txt'
    np.savetxt(file_out, data, fmt='%.6f', delimiter=' ', header='x y pdf(x,y)')
    return data


def plot_map(x, y, z):
    plt.tricontourf(x, y, z)
    plt.colorbar(label='pdf(x,y)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D PDF')
    plt.show()


def smooth_2D(values_x, values_y, nbins=1000):
    """
    Function that estimate PDF based on the sample using a gaussian kernel density estimation.
    It returns the smoothed PDF.
    """
    limits_x = np.percentile(values_x, [2, 98])
    limits_y = np.percentile(values_y, [2, 98])
    xcenters = np.linspace(limits_x[0], limits_x[1], nbins)
    ycenters = np.linspace(limits_y[0], limits_y[1], nbins)

    points = np.array([[x, y] for x in xcenters for y in ycenters])

    pdf_spread_smooth = fastkde.fastKDE.pdf_at_points(values_x, values_y, list_of_points=points)

    return pdf_spread_smooth, points[:, 0], points[:, 1]


def plot_pdf_2D(pdf_spread_smooth, x, y, parameter_2D, prefix_output):

    z = pdf_spread_smooth

    x_values = np.sort(np.unique(x))
    y_values = np.sort(np.unique(y))

    dx = x_values[1] - x_values[0]
    dy = y_values[1] - y_values[0]

    Z = z.reshape(len(y_values), len(x_values))
    X = x.reshape(len(y_values), len(x_values))
    Y = y.reshape(len(y_values), len(x_values))

    x_min = X.min()
    x_max = X.max()

    y_min = Y.min()
    y_max = Y.max()

    px = np.sum(Z, axis=0) * dy
    py = np.sum(Z, axis=1) * dx

    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], hspace=0.05, wspace=0.05)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    #  main 2D plot
    im = ax_main.imshow(Z, origin='lower', extent=(x.min(), x.max(), y.min(), y.max()), aspect='auto')

    #  colorbar
    cbar = plt.colorbar(im, ax=ax_right)
    cbar.set_label('pdf(x,y)')

    #  top projection (p(x))
    ax_top.plot(x_values, px)
    ax_top.set_ylabel('p(x)')
    ax_top.tick_params(labelbottom=False)
    ax_top.set_xlim(x_min, x_max)
    # right projection (p(y))
    ax_right.plot(py, y_values)
    ax_right.set_xlabel('p(y)')
    ax_right.tick_params(labelleft=False)
    ax_right.set_ylim(y_min, y_max)
    ax_main.set_xlabel(f'{parameter_2D[0]}')
    ax_main.set_ylabel(f'{parameter_2D[1]}')
    ax_main.set_xlim(x_min, x_max)
    ax_main.set_ylim(y_min, y_max)

    ax_top.set_title('2D PDF with 1D projections')
    file = prefix_output + 'pdf_' + parameter_2D[0] + '_' + parameter_2D[1] + '.png'
    print("Saving plot of 2D PDF to " + file)

    plt.savefig(file, bbox_inches='tight')
    plt.close()


def sample_2D(histogram, nbins, xedges, yedges):
    """
    Function to sample from a 2D histogram.
    """
    dx = np.diff(xedges)
    dy = np.diff(yedges)
    area = dx[:, None] * dy[None, :]
    probabilities = histogram * area
    probabilities /= np.sum(probabilities)

    flat_probabilities = probabilities.ravel()
    indices = np.random.choice(len(flat_probabilities), size=nbins**2, p=flat_probabilities)
    x_indices, y_indices = np.unravel_index(indices, histogram.shape)
    values_x = np.random.uniform(xedges[x_indices], xedges[x_indices + 1])
    values_y = np.random.uniform(yedges[y_indices], yedges[y_indices + 1])
    return values_x, values_y


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


def spread_2D(X, Y, counts, parameter):
    """
    Function to spread tails of the simulated distributions for better sampling.
    It fits a gaussian function to the histogram, multiplies its sigma by a scale factor,
    and adds the gaussian function to the initial histogram.
    """
    scale = 4.
    # Flatten everything (curve_fit needs 1D arrays)
    xdata = X.ravel()
    ydata = Y.ravel()
    zdata = counts.ravel()
    p0 = (
        np.max(counts),        # amplitude
        np.average(xdata, weights=zdata),       # x0
        np.average(ydata, weights=zdata),     # y0
        np.sqrt(np.average((xdata - np.average(xdata, weights=zdata))**2, weights=zdata)),        # sigma_x
        np.sqrt(np.average((ydata - np.average(ydata, weights=zdata))**2, weights=zdata)),        # sigma_y
    )

    pars, _ = curve_fit(Gauss_2d, (xdata, ydata), zdata, p0=p0)
    fit_z = Gauss_2d((xdata, ydata), pars[0]/2, pars[1], pars[2], pars[3] * scale, pars[4] * scale
                     ).reshape(counts.shape)

    print(('{:s} {:s}:  A= {:.2f} x0= {:.2f}, y0={:.2f}, sigma_x={:.2f},' +
          'sigma_y={:.2f}').format(
              parameter[0], parameter[1], pars[0], pars[1], pars[2], pars[3]*scale, pars[4]*scale))
    spread_histogram = np.add(counts, fit_z)
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


def plot_2D(xcenters, ycenters, histogram, parameter, prefix_output):
    plt.imshow(histogram, origin='lower', extent=(xcenters[0], xcenters[-1], ycenters[0], ycenters[-1]), aspect='auto')
    plt.colorbar(label='pdf(x,y)')
    plt.xlabel(parameter[0])
    plt.ylabel(parameter[1])
    file = prefix_output+parameter[0]+'_'+parameter[1]+'.png'
    print("Saving plot of used prior on " + parameter[0] + " and " + parameter[1] + " to " + file)
    plt.savefig(file)
    plt.close()


if __name__ == '__main__':

    simulated_file = 'data/OB03235/OB03235_genulens.out'
    prefix_output = 'data/OB03235/OB03235_prior_'
    data = get_data(simulated_file)
    weights = data['wtj']

    # choose parameters
    parameters = ['t_E', 'pi_E_E', 'pi_E_N']
    parameters_2D = [['pi_E_E', 'pi_E_N']]
    # convesion from MulensModel to genulens names
    conve = {
        'pi_E_E': 'pi_EE',
        'pi_E_N': 'pi_EN',
        }

    # for parameter in parameters:
    #     parameter_genulens = conve.get(parameter, parameter)
    #     values = data[parameter_genulens]
    #     pdf_orginal, pdf_spread_smooth, bins, limits = get_pdf(values, weights, parameter)
    #     plot_prior(parameter, limits, pdf_orginal, pdf_spread_smooth, prefix_output)
    #     file_out = prefix_output + parameter + '.txt'
    #     combined = np.column_stack((bins, pdf_spread_smooth(bins)))
    #     np.savetxt(file_out, combined, fmt='%.6f', delimiter=' ', header='# x pdf(x)')

    for parameter_2D in parameters_2D:
        parameter_genulens_1 = conve.get(parameter_2D[0], parameter_2D[0])
        parameter_genulens_2 = conve.get(parameter_2D[1], parameter_2D[1])
        values_x = data[parameter_genulens_1]
        values_y = data[parameter_genulens_2]
        weights = data['wtj']
        _ = get_pdf_2D(values_x, values_y, weights, parameter_2D, prefix_output)
