### To be discussed:

- periodic variables
- log10 variables etc.
- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???

### Documentation:
- Delta t_0
- binary source
- x_caustic_in etc.
- ob08092-o4_prior.yaml
- posterior files
- print model
- example with add_2450000: False

# NOW - plot trajectory:
- test if works in ulens_model_plot.py
- add satellite trajectory
- add legend
- note in README.md

## List of task to be done:

( **boldface** - do this before sending e-mail around)

- **some documentation - see above**
- Mroz+20 - finish
- MN: add option to plot best model from each mode
- MN: add more parameters to _parse_fitting_parameters_MultiNest(): n_clustering_params, max_iter, resume [previous run], const_efficiency_mode, wrapped_params [list of 0 or 1 (1 for wrap arround)], mode_tolerance, log_zero, seed [random no. generator seed], verbose [need update on sampling progress?]; FOR MORE INFO SEE: https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest/run.py AND https://github.com/farhanferoz/MultiNest/blob/master/MultiNest_v3.12/nested.F90 AND https://github.com/JohannesBuchner/MultiNest/blob/master/README
- example with add_2450000: False
- script and MM versions should be printed
- EMCEE: we should have settings in one dict - similarly to self._kwargs_MultiNest
- EMCEE backend - https://emcee.readthedocs.io/en/stable/user/backends/#emcee.backends.HDFBackend
- add one more fitting method? scipy.optimize, ultranest, https://lmfit.github.io/lmfit-py/, sfit by Jen, ???
- add check if 't_0' is covered by data and give warning if not
- print fixed parameters at begin or "no fixed parameters", so that full model can be extracted without the input file
- LD coeffs - there should be check who bands there compare to the ones in datasets
- random seed - first just print it early on (if used in calculations); then allow setting it for exact reproduction of results
- all_parameters in _get_parameters_ordered() and _check_fixed_parameters() - combine in a single one
- note that parameters are re-ordered (maybe in future add option for specifying order)
- datasets - guessing 245/246; plotting as well
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- allow plotting multiple models
- allow plotting many random models from posterior
- add beta distribution to allowed distributions (search for "gauss")
- for plot script add printing chi2 and fluxes
- EMCEE: starting parameters are read from file (e.g. 100 lines by 7 columns)
- EMCEE: some of the starting values are calculated based on equation given in yaml file, eg. "s: equation 100 / t_E" and then substitute each value of t_E and then use: "exec('out = 100 / 20.12345')" and use variable 'out'; This requires import from math of log, log10, arcsin etc.; make sure "s" in x_caustic_in is not replaced etc.; 
- if Cassan08 paramaterization is used then make sure times are >2450000.
- add automatic "obvious" checks on parameters: t_E>0, rho>0, s>0, 1>q>0 - even if they are not provided, then model should be rejected and warning given
- if magnification calculations break then give warning, reject the model, and continue
- binary source models - print fluxes of both sources separately
- warnings if plots will overwrite existing files
- check if output files (including plots) exists at the begin - similar to _warn_about_existing_output_files_MultiNest()
- plot title
- make plots tighter, i.e., reduce white space
- EMCEE: add ln_prior values to blob? At some point we will want to save that information in output files
- settings['input_file_root'] = input_file_root - in final function and use it for default output files names
- EMCEE: posterior output: 1) add log(prior), 2) add chi2 or equivalent, 3) add option to add fluxes
- EMCEE: print number of models calculated
- MN: for multimode version add option to print statistics of all modes combined
- MN: separate corner plot for each mode (requires same shift to be used in _shift_t_0_in_samples())
- periodic variables - suggest it for alpha, x_caustic_X
- check if data files exist
- allow log10() of parameter
- Event.get_chi2() - add fit_blending=False option (actually this is different in MM v2)
- allow turning off flux printing
- EMCEE: in self._sampler.run_mcmc() and option for progress=True
- warnings on time plotting and data limits - checks for add/subtract 245/246
- if code fails during fitting, then it should still print the best model found so far - add try/except in _run_fit()
- example how to run fits on a grid of (s,q)
- allow periodic (either based on number of steps, or execution time) print of best model etc.
- print every n-th model
- for parallax models check if t_0_par is fixed and give warning, if not
- fits with 0 blending flux for some datasets
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data so that anomalies are not missed etc.
- add option to adjust Y scale to plot model fully
- in _parse_fit_constraints_prior() add a check if the priors are defined for fit parameters
- flux constraints for binary source models (note that for plotting it is now set to first dataset)
- allow Model.set_magnification_methods_parameters()
- triangle and trace plots - add option to plot burn-in as well
- methods - if only single string is provided, then this is a default method
- move _get_weighted_percentile() to a separate file with utils because it doesnt depend on self; maybe there are other similar functions
- allow LD parameters to be fitted
- for trace and triangle plots, the setting `shift t_0` is common - it should be checked if it`s not set twice to different values
