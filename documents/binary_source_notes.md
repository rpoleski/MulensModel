* logically, the 2 sources will be specified in ModelParameters
* current UC: have t_0_1 and t_0_2 instead of t_0, same for u_0
* think how we calculate effective magnification - is q_f ModelParameters property? Probably No, because it can be fitted via regression

Hence, the best approach I see is to have in Model and ModelParameters private properties that will be instances of these classes but for single sources. Then the main part of Model.magnification will look something like:

```python
mag_0 = magnification_curve_0.magnification
mag_1 = magnification_curve_1.magnification
magnification = (mag_0 * f_s_0 + mag_1 * f_s_1) / (f_s_0 + f_s_1)
return magnification
```

High level functions that need changes:

* Model:
  * magnification()
  * data_magnification()
  * get_data_magnification()
  * get_ref_fluxes()
  * get_residuals()
  * functions for setting methods: set_default_magnification_method(), set_magnification_methods(), and set_magnification_methods_parameters()
  * all plotting functions
* Event:
  * get_chi2()
  * get_chi2_for_dataset()
  * get_chi2_per_point()
  * get_ref_fluxes()
  * all plotting functions (they call Model functions

