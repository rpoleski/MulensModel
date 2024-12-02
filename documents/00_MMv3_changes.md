# MulensModel V3

There were two major goals for MulensModel V3:
- Allow for N sources (However, xallarap can only be used with binary source models)
- Allow the user a direct access to partial derivatives of the magnification relative to the point lens model parameters (in addition to chi2 derivatives)

Implementing these changes involved some major updates to the internal 
architecture, but these should mostly be invisible to the user (although they 
add new functionality).

We have also cleaned up some problems/bugs/features that will impact the user experience:
- Shifted the convention of `alpha` by 180 deg to make MM compatible with convention of Skowron et al. 2011.
- Removing `astropy` units from microlensing parameters.
- Generally improving consistency in naming conventions.
- Deprecated functions and properties removed.

A detailed list of changes follows.

## Changes to Major Classes

### Event Class

- REMOVED deprecated functions and attributes:
  - `get_ref_fluxes()`: REMOVED keywords `data_ref` and `fit_blending`.
  - REMOVED keyword `fit_blending` of methods `get_ref_fluxes()`, `get_chi2()`, 
`get_chi2_for_dataset()`, `get_chi2_per_point()`, and `get_chi2_gradient()`. 
Use `fix_blend_flux` instead.
  - `reset_best_chi2()`
  - `best_chi2`
  - `best_chi2_parameters`

### Model Class

- ADDED `get_magnification_curve()` and `get_magnification_curves()` (for multi-source models).
- `get_lc()`: updated documentation to reflect that this returns MAGNITUDES (not magnifications)
- `get_magnification()`: 
   - REMOVED keyword `flux_ratio_constraint`
   - default value of `separate` changed to None (which defaults to `True` if `source_flux_ratio` is provided and to `False` otherwise).
- REMOVED deprecated functions and attributes:
  - `plot_magnification()`: REMOVED keyword `flux_ratio_constraint`.
  - `plot_lc()`: REMOVED keywords `data_ref`, `flux_ratio_constraint`, `fit_blending`, `f_source`, and `f_blend`.
  - `plot_trajectory()`: REMOVED keyword `show_data`.
  - `set_default_magnification_method()`. Use `Model.default_magnification_method = XXX`, instead.
  - `magnification()` replaced by `get_magnification()`.
  - `reset_plot_properties()`

### ModelParameters Class

Even if the user was initializing `Model` with a `dict`, that `dict` is 
immediately converted into a `ModelParameters` object, so many of these changes 
could affect the user experience, even if they had not been interacting directly 
with the `ModelParameters` class.

- `Astropy` units REMOVED (i.e., now everything is a `float`) from properties: `t_star`, `t_eff`, `t_E`, `alpha`, `ds_dt`, `dalpha_dt`, `t_star_1`, `t_star_2`, `gamma_parallel`, `gamma_perp`, and `gamma`. Also REMOVED from function `get_alpha()`.

- Property `pi_E` REMOVED. Properties `pi_E_E` and `pi_E_N` are still there.
- Properties that are not defined now raise `AttributeError` (previously `rho` was returning None, and some other ones raised `KeyError`).
- Warnings for unexpected values of angles (e.g., 1000 deg) are not raised anymore.
- Changes in default values for `t_0_par` and `t_0_kep`:
  - If `t_0_par` is not defined, then it defaults to `t_0_kep`, `t_0`, `t_0_1` in that order. 
  - If `t_0_kep` is not defined, then it defaults to `t_0_par`, `t_0`, `t_0_1` in that order.

Also, REMOVED `which_parameters()` from modelparameters.py because it not very 
useful and it was very complicated.

### MulensData Class

- REMOVED deprecated functions and attributes:
  - coordinates (i.e., `coords`, `ra`, `dec`) are no longer associated with 
`MulensData` objects.
  - REMOVED `model` and `plot_residuals` keywords from `plot()`.
  
### Trajectory Class

ADDED new keyword arguments to `__init__`: `x` and `y`. They can be provided 
instead of the `times` argument.

## Changes to Other Classes

### FitData Class

- Now supports `fit_fluxes` for multiple sources.
- ADDED:
  - `get_d_A_d_rho()`
  - `magnification_curve` and `magnification_curves`
- DEPRECATED subclass `FSPL_Derivatives`. Most its methods are replaced by `PointLens` class(es)

- DEPRECATED:
  - `get_d_A_d_u_for_PSPL_model()` 
  - `get_d_A_d_u_for_FSPL_model()` 
  - `get_d_A_d_u_for_point_lens_model()`

All `get_d_A_d_u*` and `get_d_u_d_params*` functions got moved to the 
`PoinLens*Magnification` classes in pointlens.py (see below). If desired, they 
could be accessed instead through `FitData.magnification_curve.XXX()`.

- REMOVED deprecated functions and attributes:
  - In `get_residuals`, REMOVED `type` keyword. Use `phot_fmt` instead.

### Point Lens
The PointLens and BinaryLens classes were subdivided into separate classes for 
each magnification method.

- Function `get_pspl_magnification()` from `pointlens.py` was removed. Use classes listed in next point instead.
- Point Lens classes ADDED:
  - `PointSourcePointLensMagnification` 
  - `FiniteSourceUniformGould94Magnification` 
  - `FiniteSourceLDYoo04Magnification` 
  - `FiniteSourceUniformWittMao94Magnification` 
  - `FiniteSourceLDWittMao94Magnification` 
  - `FiniteSourceUniformLee09Magnification`
  - `FiniteSourceLDLee09Magnification`
- `PointLensWithShear` REPLACED by `PointSourcePointLensWithShearMagnification`.

### Binary Lenses
- `BinaryLens` class REMOVED.
- Binary classes ADDED: 
  - `BinaryLensPointSourceWM95Magnification` 
  - `BinaryLensPointSourceVBBLMagnification` 
  - `BinaryLensPointSourceMagnification` 
  - `BinaryLensQuadrupoleMagnification` 
  - `BinaryLensHexadecapoleMagnification` 
  - `BinaryLensVBBLMagnification` 
  - `BinaryLensAdaptiveContouringMagnification` 
  - `BinaryLensPointSourceWithShearWM95Magnification` 
  - `BinaryLensPointSourceWithShearVBBLMagnification`
- `Caustics` REPLACED by `CausticsBinary`.
- `CausticsWithShear` REPLACED by `CausticsBinaryWithShear`.

### Other
- `PointLensFiniteSource` REPLACED by `B0B1Utils`. 
- ADDED `EllipUtils`.
- `MagnificationCurve`:
  - ADDED methods `get_d_A_d_params()` and `get_d_A_d_rho()`.
  - ADDED property `methods_indices`.
