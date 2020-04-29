# High-level script for fitting microlensing models with MulensModel

This script aims at allowing easy access to fitting capabilities that gives MulensModel.
All settings are passed via YAML files, which are human- and machine-readable. 

Example usage:

```
python ulens_model_fit.py ob08092-o4_minimal.yaml
```

should produce fitted model parameters in a few seconds. Please have a look at `ob08092-o4_minimal.yaml` - it has only basic settings. In many cases, one can fit a reasonable point-source point-lens model by jus changing file name and mean value of `t_0`.

More complicated example that also produced plot of the best model with residuals and the triangle plot:

```
python ulens_model_fit.py ob08092-o4.yaml
```

There are many more features to be added - _please let authors know what specific needs you have_.

Please note that this code is a high-level example for MulensModel, but it uses fitting algorithms that are not part of MulensModel. The latter allows many microlensing calculations including chi^2 for given data and model parameters, but does not have built-in fitting capabilities.

