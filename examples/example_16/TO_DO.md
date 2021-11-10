### To be discussed:

- periodic variables
- log10 variables etc.
- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???
- how the output should look like?

### Documentation:
- Delta t_0
- binary source
- x_caustic_in etc.
- ob08092-o4_prior.yaml
- posterior files
- print model

## List of task to be done:

( **boldface** - do this before sending e-mail around)

- **some documentation - see above**
- **add one more fitting method?**
- **requirements.txt** - corner >= 2.0.0
- Mroz+20 - finish
- trace plot (give option to include burn-in)
- script and MM versions should be printed
- EMCEE backend - https://emcee.readthedocs.io/en/stable/user/backends/#emcee.backends.HDFBackend
- add check if 't_0' is covered by data and give warning if not
- print fixed parameters at begin or "no fixed parameters", so that full model can be extracted without the input file
- LD coeffs - there should be check who bands there compare to the ones in datasets
- random seed - first just print it early on (if used in calculations); then allow setting it for exact reproduction of results
- all_parameters in _get_parameters_ordered() and _check_fixed_parameters() - combine in a single one
- note that parameters are re-ordered (maybe in future add option for specifying order)
- datasets - guessing 245/246; plotting as well
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- allow plotting multiple models
- add beta distribution to allowed distributions (search for "gauss")
- for plot script add printing chi2 and fluxes
- starting parameters are read from file (e.g. 100 lines by 7 columns)
- some of the starting values are calculated based on equation given in yaml file, eg. "s: equation 100 / t_E" and then substitute each value of t_E and then use: "exec('out = 100 / 20.12345')" and use variable 'out'; This requires import from math of log, log10, arcsin etc.; make sure "s" in x_caustic_in is not replaced etc.; 
- if Cassan08 paramaterization is used then make sure times are >2450000.
- add automatic "obvious" checks on parameters: t_E>0, rho>0, s>0, 1>q>0 - even if they are not provided, then model should be rejected and warning given
- if magnification calculations break then give warning, reject the model, and continue
- binary source models - print fluxes of both sources separately
- Fitting method to be added: scipy.optimize, pymultinest, ultranest, https://lmfit.github.io/lmfit-py/, ???
- allow plotting many random models from posterior
- warnings if plots will overwrite existing files
- check if output files (including plots) exists at the begin
- plot title
- make plots tighter, i.e., reduce white space
- Add ln_prior values to blob? At some point we will want to save that information in output files
- settings['input_file_root'] = input_file_root - in final function and use it for default output files names
- posterior output: 1) add log(prior), 2) add chi2 or equivalent, 3) add option to add fluxes
- print number of models calculated
- periodic variables - suggest it for alpha, x_caustic_X
- check if data files exist
- allow log10() of parameter
- Event.get_chi2() - add fit_blending=False option (actually this is different in MM v2)
- allow turning off flux printing
- warnings on time plotting and data limits - checks for add/subtract 245/246
- if code fails during fitting, then it should still print the best model found so far - add try/except in _run_fit()
- example how to run fits on a grid of (s,q)
- allow periodic print of best model etc.
- print every n-th model
- plot trajectory
- for parallax models check if t_0_par is fixed and give warning, if not
- fits with 0 blending flux for some datasets
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data so that anomalies are not missed etc.
- add scipy to _check_imports() - requires siginificant code to be added to _check_imports() in order to find out if t_E prior is used
- in _parse_fit_constraints_prior() add a check if the priors are defined for fit parameters
- flux constraints for binary source models (note that for plotting it is now set to first dataset)
- allow Model.set_magnification_methods_parameters()
- methods - if only single string is provided, then this is a default method
- print current best model - each minute, each nth model etc.
- allow LD parameters to be fitted
