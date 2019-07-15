"""
Prepare files to be used as tests.
"""
import numpy as np

from MulensModel import Model


MAG_ZEROPOINT = 18.


def simulate_PSPL(file_name, n_data, t_start=None, t_stop=None, u_0=None,
                  magnification_function=None):
    """simulate PSPL light curve and save to file"""
    t_0 = 2456900.
    t_E = 20.
    u_0 = 0.01
    if t_start is None:
        t_start = t_0 - 4. * t_E
    if t_stop is None:
        t_stop = t_0 + 4. * t_E
    relative_sigma = 0.03
    add_magnitude_error = 0.003
    flux_source = 3.
    flux_blend = 0.1

    times = np.sort(np.random.uniform(t_start, t_stop, n_data))
    if magnification_function is None:
        tau = (times - t_0) / t_E
        u2 = tau**2 + u_0**2
        magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    else:
        magnification = magnification_function(times)
    flux = magnification * flux_source + flux_blend
    diff = np.random.normal(0., 1., n_data)
    sigma = flux * relative_sigma
    flux_observed = flux + sigma * diff
    mag_observed = MAG_ZEROPOINT - 2.5 * np.log10(flux_observed)
    sigma_observed = np.sqrt((sigma/flux_observed)**2 + add_magnitude_error**2)

    data_out = np.array([times, mag_observed, sigma_observed]).T
    # data_out = np.array([times, magnification, sigma_observed]).T
    np.savetxt(file_name, data_out)


if __name__ == '__main__':
    simulate_PSPL('test_100.txt', 100)
    simulate_PSPL('test_1000.txt', 1000)
    simulate_PSPL('test_10000.txt', 10000)

    parameters = {'t_0': 2456900., 'u_0': 0.1, 't_E': 50.,
                  'pi_E_N': 0.6, 'pi_E_E': 0.8}
    model = Model(parameters, coords="18:00:00.00 -30:00:00.0")
    model.parallax(earth_orbital=True)

    kwargs = {'magnification_function': model.magnification,
              't_start': parameters['t_0']-80.,
              't_stop': parameters['t_0']+80.}

    simulate_PSPL('test_100_piE.txt', 100, **kwargs)
    simulate_PSPL('test_1000_piE.txt', 1000, **kwargs)
    simulate_PSPL('test_10000_piE.txt', 10000, **kwargs)
