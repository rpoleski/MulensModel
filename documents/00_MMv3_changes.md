# MulensModel v3

## Major changes from the user perspective:

- Changing convention of alpha angle, i.e., shifting by 180 deg. This makes MM compatible with convention of Skowron et al. 2011.
- Allowing more than 2 sources.
- Removing astropy units from microlensing parameters.
- Deprecated functions removed
- something about derivatives

## Additional Changes

### New/Renamed Classes
- PointLensFiniteSource replaced by B0B1Utils. Introduction of EllipUtils.
- BinaryLens class removed.
- Binary classes introduced: BinaryLensPointSourceWM95Magnification, BinaryLensPointSourceVBBLMagnification, BinaryLensPointSourceMagnification, BinaryLensQuadrupoleMagnification, BinaryLensHexadecapoleMagnification, BinaryLensVBBLMagnification, BinaryLensAdaptiveContouringMagnification, BinaryLensPointSourceWithShearWM95Magnification, and BinaryLensPointSourceWithShearVBBLMagnification.
- Caustics replaced by CausticsBinary.
- CausticsWithShear replaced by CausticsBinaryWithShear.

### Technical Changes to Major Classes

#### Event Class

- Event.get_ref_fluxes() does not take arguments. 

#### Model Class

#### MulensData Class

### Additional Changes


The parameter fit_blending of Event methods get_chi2(), get_chi2_for_dataset(), get_chi2_per_point(), and get_chi2_gradient() removed. Use fix_blend_flux instead.
Parameter type removed from FitData.get_residuals().
Methods of class FitData removed: get_d_A_d_u_for_PSPL_model(), get_d_A_d_u_for_FSPL_model(), and get_d_A_d_u_for_point_lens_model().
Properties FitData.magnification_curve and magnification_curves removed.
Class FSPL_Derivatives removed. Most its methods are replaced by PointLens class.
Methods MagnificationCurve.get_d_A_d_params() and get_d_A_d_rho() added.
Property MagnificationCurve.methods_indices added.
Parameter flux_ratio_constraint of Model.plot_magnification() removed.
Parameters of data_ref, flux_ratio_constraint, fit_blending, f_source, and f_blend of Model.plot_lc() removed.
Parameter show_data of Model.plot_trajectory() removed.
Model.get_magnification(): parameter flux_ratio_constraint remove and parameter separate default value changed to None (which defaults to True if source_flux_ratio is provided and to False otherwise).
Model.magnification() replaced by Model.get_magnification().
Model.get_magnification_curve() and get_magnification_curves() added.
Function which_parameters() from modelparameters.py removed as not useful and very complicated.
Property ModelParameters.pi_E removed. Properties pi_E_E and pi_E_N are still there.
Astropy units removed (i.e., now everything is a float) from ModelParameters properties: t_star, t_eff, t_E, alpha, ds_dt, dalpha_dt, t_star_1, t_star_2, gamma_parallel, gamma_perp, and gamma. Also removed from function get_alpha().
In ModelParameters, the properties that are not defined raise AttributeError (previously rho was returning None, and some other ones raised KeyError).
Warnings for unexpected values of angles (e.g., 1000 deg) in ModelParameters are not raised anymore.
Changes in default values for t_0_par and t_0_kep. If the t_0_par is not defined, then it defaults to t_0_kep, t_0, t_0_1 in that order. If the t_0_kep is not defined, then it defaults to t_0_par, t_0, t_0_1 in that order.
MulensData arguments coords, ra, and dec removed. Also property coords was removed.
Parameters model and plot_residuals removed from MulensData.plot().
Function get_pspl_magnification() from pointlens.py was removed. Use classes listed in next point instead.
Added classes PointSourcePointLensMagnification, FiniteSourceUniformGould94Magnification, FiniteSourceLDYoo04Magnification, FiniteSourceUniformWittMao94Magnification, FiniteSourceLDWittMao94Magnification, FiniteSourceUniformLee09Magnification, and FiniteSourceLDLee09Magnification.
PointLensWithShear replaced by PointSourcePointLensWithShearMagnification.
Added new arguments of Trajectory: x and y. They can be provided instead of the times argumen