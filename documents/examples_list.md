Tutorials MulensModel:

* [Basic usage tutorial](https://rpoleski.github.io/MulensModel/tutorial.html),
* [Fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html),
* [Microlensing parallax fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pi_E.html),
* [**Explanation of microlensing parameters**](parameter_names.pdf),
* [**Explanation of methods used for calculating magnification**](magnification_methods.pdf),
* [Instructions on getting satellite positions](Horizons_manual.md) - useful only if you have satellite data.

[Examples on how to use the code](../examples/):
* [Example 01](../examples/example_01_models.py) - plot simple point-source/point-lens (PSPL) model and model with planetary lens,
* [Example 02](../examples/example_02_fitting.py) - fit PSPL model to the data using [scipy.optimize.minimize()](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html),
* [Example 03](../examples/example_03_mulenssystem.py) - define PSPL model using physical properties and plot the resulting magnification curve,
* [Example 04](../examples/example_04_einsteinring.py) - calculate the Einstein ring size for a grid of lens masses and distances,
* [Example 05](../examples/example_05_MB08310.py) - plot multiple datasets for a single model, plot the residuals, and do this both in magnitude and magnification spaces,
* [Example 06](../examples/example_06_fit_parallax_EMCEE.py) - fit parallax model using [EMCEE](https://emcee.readthedocs.io/en/stable/),
* [Example 07](../examples/example_07_fit_parallax_MN.py) - fit parallax model using [MultiNest](https://github.com/JohannesBuchner/PyMultiNest),
* [Example 08](../examples/example_08_planet_grid_fitting.ipynb) - shows how to fit simulated WFIRST light curve with planetary model,
* [Example 09](../examples/example_09_gradient_fitting.py) - fit point lens model using chi^2 gradient,
* [Example 10](../examples/example_10_fitting_and_fluxes.py) - fit model and extract posterior fluxes, use [config file](../examples/example_10.cfg) to pass all parameters,
* [Example 11](../examples/example_11_binary_source.py) - simulate and fit binary source event,
* [Example 12](../examples/example_12_fit_satellite_parallax_EMCEE.py) - fit parallax model to ground and satellite data, plot models and trajectories at the end,
* [Example 13](../examples/example_13_caustic_sampling.py) - fit planetary event using caustic entrance and exit epochs as parameters (uses [config file](../examples/example_13.cfg)),
* [Example 14](../examples/example_14_caustic_plotting.py) - plot caustic using standard method and uniform sampling,
* [Example 15](../examples/example_15_fitting.py) - fitting binary lens model with many options - use config file for [ob05390](../examples/example_15_ob05390_v1.cfg) or [mb07192](../examples/example_15_mb07192_v1.cfg); settings are read by [this file](../examples/example_15_read.py),
* [Example 16](../examples/example_16/) - **high-level fitting example** where all settings are read from a human-readable YAML file, there is a separate description of that example in [this README file](../examples/example_16/README.md),
* [Example 17](../examples/example_17_1L2S_plotting.py) - plotting binary source model,
* [Example 18](../examples/example_18_simulate.py) - simulate a light curve and save it; example input files: [example_18_input_1.yaml](../examples/example_18_input_1.yaml) and [example_18_input_2.yaml](../examples/example_18_input_2.yaml),
* [Example 19](../examples/19_binary_source_binary_lens.py) - make plots of a model that has 2 sources and 2 lenses,
* Three files producing plots presented in paper describing MulensModel: [plots_1.py](../examples/plots_1.py), [plots_2.py](../examples/plots_2.py), and [plots_3.py](../examples/plots_3.py).

[MulensModel documentation](https://rpoleski.github.io/MulensModel/) includes description of input and output of every function. 

If you have used MulensModel and wrote a code that you think might be useful for others, then please send it to code authors or submit a pull request.

