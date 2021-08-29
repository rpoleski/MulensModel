"""
Script for simulating microlensing lightcurves.
All settings are controlled from a yaml file. Example input:

    python example_18_simulate.py example_18_input_1.yaml

    python example_18_simulate.py example_18_input_2.yaml

The first one is very simple. The second one is more complicated.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import yaml

import MulensModel as mm


def simulate_lc(
        parameters, time_settings, file_out,
        coords=None, methods=None,
        flux_source=1000., flux_blending=0.,
        relative_uncertainty=0.01,
        plot=True, subtract_2450000=True):
    """
    Simulate and save light curve.

    Parameters :
        parameters: *dict*
            Parameters of the model - keys are in MulensModel format, e.g.,
            't_0', 'u_0', 't_E' etc.

        time_settings: *dict*
            Sets properties of time vector. It requires key `type`, which can
            have one of two values:
            - `random` (requires `n_epochs`, `t_start`, and `t_stop`) or
            - `evenly spaced` (settings passed to `Model.set_times()`).

        file_out: *str*
            Name of the file to be saved.

        coords: *str*
            Event coordinates for parallax calculations, e.g.,
            "17:34:51.15 -30:29:28.27".

        methods: *list*
            Define methods used to calculate magnification. The format is
            the same as MulensModel.Model.set_magnification_methods().

        flux_source: *float*
            Flux of source.

        flux_blending: *float*
            Blending flux.

        relative_uncertainty: *float*
            Relative uncertainty of the simulated data (this is close to
            sigma in magnitudes).

        plot: *bool*
            Plot the data and model at the end?

        subtract_2450000: *bool*
            Do you want shorter JD values?
    """
    model = mm.Model(parameters, coords=coords)

    if time_settings['type'] == 'random':
        raw = np.random.rand(time_settings['n_epochs'])
        dt = time_settings['t_stop'] - time_settings['t_start']
        times = time_settings['t_start'] + np.sort(raw) * dt
    elif time_settings['type'] == 'evenly spaced':
        times = model.set_times(**time_settings)
    else:
        raise ValueError("unrecognized time_settings['type']: " +
                         time_settings['type'])

    if methods is not None:
        model.set_magnification_methods(methods)

    magnification = model.get_magnification(times)

    flux = flux_source * magnification + flux_blending
    flux_err = relative_uncertainty * flux

    flux *= 1 + relative_uncertainty * np.random.randn(len(flux))

    data = mm.MulensData([times, flux, flux_err], phot_fmt='flux')
    event = mm.Event([data], model)
    print("chi^2: {:.2f}".format(event.get_chi2()))

    if subtract_2450000:
        subtract = 2450000.
    else:
        subtract = 0.

    np.savetxt(file_out,
               np.array([times-subtract, data.mag, data.err_mag]).T,
               fmt='%.4f')

    if plot:
        model.plot_lc(subtract_2450000=subtract_2450000,
                      t_start=np.min(times), t_stop=np.max(times),
                      source_flux=flux_source, blend_flux=flux_blending)
        data.plot(phot_fmt='mag', subtract_2450000=subtract_2450000)
        plt.show()


if __name__ == '__main__':
    input_file = sys.argv[1]
    with open(input_file, 'r') as in_file:
        settings = yaml.safe_load(in_file)

    simulate_lc(**settings)
