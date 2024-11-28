# High-level script for fitting microlensing models with MulensModel

First things first:

__This script aims at allowing easy access to fitting microlensing models with the [MulensModel package](https://github.com/rpoleski/MulensModel) and it's not the best place to learn the [MulensModel](https://github.com/rpoleski/MulensModel) itself.__

Allowing easy access to many functions results in somehow complicated code, so if you want to learn MulensModel usage, then we suggest to start with other examples. Also please note that MulensModel has more capabilities than provided in this script.

Here all settings are passed via YAML files, which are human- and machine-readable. If the script syntax is unclear, then please search for information on YAML format files (or see [this link](https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html)).

### Install required packages

I suggest to start with:
```
pip install -r requirements.txt
```
so that you have all required packages.

### Basic usage

In this and following few sections I show how to fit model using [EMCEE](https://emcee.readthedocs.io/en/stable/) implementation of [MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo). The other possibility is to use [MultiNest](http://johannesbuchner.github.io/PyMultiNest/) and it's presented at the end.

Example usage:

```python
python ulens_model_fit.py ob08092-o4_minimal.yaml
```

should produce fitted model parameters in a few seconds. Please have a look at [`ob08092-o4_minimal.yaml`](ob08092-o4_minimal.yaml) - it has only basic settings. In many cases, one can fit a reasonable point-source point-lens model by just changing file name and mean value of `t_0`.

This minimal script does not plot the model. You can plot the model at the end of fitting, or run a separate script ([`ulens_model_plot.py`](ulens_model_plot.py)) that only makes the plot. You have to take 3rd and 4th last line from output of the above script, add them in `model` secion of the YAML file and add the information where you want your plot to be saved:

```python
python ulens_model_plot.py ob08092-o4_minimal_plot.yaml
```

This should produce file `ob08092-o4_minimal_plot.png`. Compare the two above YAML files to see the differences. You can remove three sections from plotting YAML file, because they're ignored anyway.

More complicated example that will also produce plots of the best model with residuals and the triangle plot:

```python
python ulens_model_fit.py ob08092-o4.yaml
```

In [`ob08092-o4.yaml`](ob08092-o4.yaml) you can see how format of these YAML files mirrors MulensModel API - see the second line and compare it to [MulensData API](https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html).

### Binary lens fitting

You can specify the methods used for calculating magnification. For example, fit the first microlensing planet (calculations may take a few minutes):

```python
python ulens_model_fit.py ob03235_1.yaml
```

Note that this YAML file will result in a warning message. The message is caused by the fact that in some cases, flux minus its uncertainty results in negative values which cannot be translated to magnitudes when plotting. The warning is related to only plotting, not fitting. You can ignore this warning.

There are many more features to be added - _please let authors know what specific needs you have_.

Please note that this code is a high-level example for MulensModel, but it uses fitting algorithms that are not part of MulensModel. The latter allows many microlensing calculations including chi^2 for given data and model parameters, but does not have built-in fitting capabilities.


### Annual parallax fitting

Let's try to fit the parallax model for OB05086. First, fit model without parallax:

```python
python ulens_model_fit.py ob05086_1.yaml
```

We see that some of the points are not well fit. Hence, we will try to fit the parallax model. First we update starting points based on results from the first fit (e.g., `t_E = 100`). We also add parallax parameters: `pi_E_N` and `pi_E_E`. We limit both of them to `(-1, 1)` range. There is one more piece of information definitely needed: event coordinates. These are: `18:04:45.71 -26:59:15.2`. In YAML file it's under `model` and `coords`. Finally, we want to set parameter reference epoch: `t_0_par`. The fit will be much much slower without it. We choose an epoch close to `t_0`. This is not a fitting parameter, so it's placed in `fixed_parameters`. All these changes are in `ob05086_2.yaml` file. Please compare it previous file to see all changes. And run:

```python
python ulens_model_fit.py ob05086_2.yaml
```

We see that chi^2 improved significantly - by more than 400. This clearly indicates the parallax model is better. We see it also on the lightcurve plot. However there is a problem. The blending flux is significantly negative, which is unphysical. Let's try degenerate solution, i.e., we change `u_0` sign:

```python
python ulens_model_fit.py ob05086_3.yaml
```

The chi^2 improved very slightly over positive `u_0` parallax model, only by 2.3. In this model the blending flux is positive, so it's the best model. Congratulations!


### Priors and constraints

It is possible to specify additional fit constraints in the input file. Currently, empirical `t_E` distribution and constraining negative blending flux are allowed. To see how it works, let's go back to `ob08092-o4.yaml`. If you carefully look at the output printed, you will see that the blending flux (`flux_b_1` in output) is negative within 1-sigma. This is somehow unphysical. To prevent it, we can disfavor negative blending flux models - see [`ob08092-o4_prior_1.yaml`](ob08092-o4_prior_1.yaml):

```python
python ulens_model_fit.py ob08092-o4_prior_1.yaml
```

You will see that the blending flux changed - the median increased and uncertainties decreased. Also, positive and negative uncertainties changed to asymmetric. You can have a look at triangle plots from two runs to see the difference.

We can also specify prior on `t_E`. In [`ob08092-o4_prior_2.yaml`](ob08092-o4_prior_2.yaml) we use the results from [Mroz et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017Natur.548..183M/abstract):

```python
python ulens_model_fit.py ob08092-o4_prior_2.yaml
```

In this case, the parameters don't change much because `t_E` is well constrained by the data.

For detailed description of different options, see `fit_constraints` in [ulens\_model\_fit.py](ulens_model_fit.py).


### Fit using pyMultiNest

The pyMultiNest is one of the implementations of [nested sampling](https://en.wikipedia.org/wiki/Nested_sampling_algorithm) - a method that has similar goals to frequently used MCMC approach. There are complicated posteriors in which EMCEE fails and pyMultiNest works without a problem. The latter method is also capable of automatically finding separate posterior modes and exploring each one separately.

The most basic usage of pyMultiNest is presented in [ob08092-o4\_minimal\_MN.yaml](ob08092-o4_minimal_MN.yaml). Note that instead of `starting_parameters` there are `prior_limits` (based on these settings it's decided which method will be used). 

More advanced pyMultiNest input file is [ob08092-o4\_MN.yaml](ob08092-o4_MN.yaml). It illustrates trivial degeneracy u0 vs. -u0. Note that `multimodal` option is turned on, so each mode is reported separately. 

### More options

There are many options and more are being added. The file [ob03235\_2\_full.yaml](ob03235_2_full.yaml) presents all options currently available:

```python
python ulens_model_fit.py ob03235_2_full.yaml
```

### Your own parametersization?

Sometimes one wants to fit using different parameters then the ones defined in MulensModel. As an example, you may be fitting a wide-orbit planet model with two separate peaks. In that case, the planet peak can be read from the light-curve easily, but it's not a standard parameter. In that case, it's enough that you define a function that translates parameters and add a few lines of code. Here is an example:

```python
python reparametrization.py reparametrization_ob08092_O3.yaml
```

Note that all other features of `ulens_model_fit.py` are available.


### More information

* Some more information on API can be found at the top of [ulens\_model\_fit.py file](https://github.com/rpoleski/MulensModel/blob/master/examples/example_16/ulens_model_fit.py) - see docstrings for UlensModelFit class. Please keep in mind that all keywords are read from YAML type file.
* Julian Dates are long numbers and which may cause problems (e.g., too many numbers to be displayed properly on axis label), hence, in this example we add 2450000 to all input data and subtract it from plots. Note that all calculations are carried out using full HJD, so e.g., Earth's positions are calculated for proper epochs.
* If you want to plot to screen then do not provide the name of output file for plot, e.g., you can remove last line in [ob08092-o4\_minimal\_plot.yaml](ob08092-o4_minimal_plot.yaml).
* In output, "Best model" is the one with the highest probability, which if different from the smallest chi2 model if informative priors are applied.
* I have many plans to add more options and capabilities. Please let me know, what you need and I'll try to make it my priority.

