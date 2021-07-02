"""
Script for plotting the model using UlensModelFit class.
All the settings are read from a YAML file.
"""
import sys
from os import path
import yaml

from ulens_model_fit import UlensModelFit


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')

    input_file = sys.argv[1]
    input_file_root = path.splitext(input_file)[0]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    # Remove settings that are not used for plotting:
    keys = ["starting_parameters", "min_values", "max_values",
            "fitting_parameters"]
    for key in keys:
        settings[key] = None
    if "plots" in settings:
        if "triangle" in settings["plots"]:
            settings["plots"].pop("triangle")

    ulens_model_fit = UlensModelFit(**settings)

    ulens_model_fit.plot_best_model()
