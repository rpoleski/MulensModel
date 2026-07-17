"""
Plot trajectory based on settings in yaml file.
"""
import matplotlib.pyplot as plt
import yaml
import sys

import MulensModel as mm


def plot_trajectory(parameters, values, time_range=None, file_out=None):
    """
    Plot trajectory.

    Parameters :
        parameters: *list* or *str*
            Parameters to be set, e.g., ['t_0', 'u_0', 't_E', 's', 'q', 'alpha']
        values: *list* or *str*
            Values of these parameters e.g., '2456789.0 0.1 100 3.14 0.5 160'
        time_range: *list* of 2 *floats* or *None*
            Time range to be plotted; set in automated way for *None*
        file_out: *str* or *None*
            Name of output file. *None* displays output on screen.
    """
    if isinstance(parameters, str):
        parameters = parameters.split()

    if isinstance(values, str):
        values = [float(v) for v in values.split()]

    model = mm.Model(dict(zip(parameters, values)))

    if time_range is None:
        t_0 = model.parameters.t_0
        t_E = model.parameters.t_E
        time_range = [t_0-1.5*t_E, t_0+1.5*t_E]

    model.plot_trajectory(t_range=time_range, caustics=True)

    if file_out is None:
        plt.show()
    else:
        plt.savefig(file_out)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('One parameter needed - yaml file')

    with open(sys.argv[1]) as in_file:
        settings = yaml.safe_load(in_file)

    plot_trajectory(**settings)
