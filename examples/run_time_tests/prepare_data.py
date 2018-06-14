"""
Prepare files to be used as tests.
"""
import numpy as np


MAG_ZEROPOINT = 18.

def simulate_PSPL(file_name, n_data):
    """simulate PSPL light curve and save to file"""
    t_0 = 2456000.
    t_E = 20.
    u_0 = 0.01
    t_start = t_0 - 4. * t_E
    t_stop = t_0 + 4. * t_E
    relative_sigma = 0.03
    add_magnitude_error = 0.003
    flux_source = 3.
    flux_blend = 0.1

    times = np.random.uniform(t_start, t_stop, n_data)
    tau = (times - t_0) / t_E
    u2 = tau**2 + u_0**2
    magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    flux = magnification * flux_source + flux_blend
    diff = np.random.normal(0., 1., n_data)
    sigma = flux * relative_sigma
    flux_observed = flux + sigma * diff
    mag_observed = MAG_ZEROPOINT - 2.5 * np.log10(flux_observed)
    sigma_observed = np.sqrt((sigma/flux_observed)**2 + add_magnitude_error**2)

    data_out = np.array([times, mag_observed, sigma_observed]).T
    np.savetxt(file_name, data_out)


if __name__ == '__main__':
    simulate_PSPL('test_100.txt', 100)
    simulate_PSPL('test_1000.txt', 1000)
    simulate_PSPL('test_10000.txt', 10000)

