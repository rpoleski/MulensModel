# High-level script for fitting microlensing models with MulensModel

First things first:

__This script aims allowing easy access to fitting microlensing models with MulensModel and it's not the best place to learn MulensModel itself.__

Allowing easy access to many functions results in somehow complicated code, so if you want to learn MulensModel usage, then we suggest to start with other examples. Also please note that MulensModel has more capabilities than provided in this script.

Here all settings are passed via YAML files, which are human- and machine-readable. If the script syntax is unclear, then please search for information on YAML format files.

### Basic usage

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

