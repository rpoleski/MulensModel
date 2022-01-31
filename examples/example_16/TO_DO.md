# test:
python ulens_model_fit.py ob08092-o4_minimal_MN.yaml

../example_07_fit_parallax_MN.py
https://github.com/JohannesBuchner/PyMultiNest
https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest/solve.py
https://github.com/JohannesBuchner/PyMultiNest/blob/master/pymultinest/run.py

# TO DO:
 - flux info - pymultinest: 1) pass Analyzer, 2) can we calculated fluxes in LogLikelihood, not Prior, 3) code changes needed
 - .gitignore
 - print info on different modes, including posterior mode probability
 - make sure different prior settings are consistent
 - _parse_fitting_parameters_MN() - we need more parameters there
 - self._return_fluxes = False is currnetly used
 - XXX
 - self._flat_priors
 - test fixed_parameters - in input file
 - give warning if files outputfiles_basename* exist - early on!
 - if there is first working version - let interested people know about this branch and ask them for input
 - requirements.txt - see below; also note it in README
 - documentation of \_\_init\_\_
 - pycodestyle
 - example full file for MN
 - add it to README, or a task to do it later
 - maybe add task: user directly says which method to fit
 - clean this task

# DONE:
 - _ln_prior
 - min/max_values cannot be set
 - output root file not provided - what to do then


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
- **add one more fitting method?** scipy.optimize, pymultinest, ultranest, https://lmfit.github.io/lmfit-py/, ???
- **requirements.txt** - corner >= 2.0.0, MM>2.0, yaml, pymultinest, ???
- Mroz+20 - finish
- script and MM versions should be printed
- EMCEE backend - https://emcee.readthedocs.io/en/stable/user/backends/#emcee.backends.HDFBackend
- plot trajectory
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
- starting parameters are read from file (e.g. 100 lines by 7 columns)
- some of the starting values are calculated based on equation given in yaml file, eg. "s: equation 100 / t_E" and then substitute each value of t_E and then use: "exec('out = 100 / 20.12345')" and use variable 'out'; This requires import from math of log, log10, arcsin etc.; make sure "s" in x_caustic_in is not replaced etc.; 
- if Cassan08 paramaterization is used then make sure times are >2450000.
- add automatic "obvious" checks on parameters: t_E>0, rho>0, s>0, 1>q>0 - even if they are not provided, then model should be rejected and warning given
- if magnification calculations break then give warning, reject the model, and continue
- binary source models - print fluxes of both sources separately
- warnings if plots will overwrite existing files
- check if output files (including plots) exists at the begin
- plot title
- make plots tighter, i.e., reduce white space
- print autocorrelation time
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
- allow periodic (either based on number of steps, or execution time) print of best model etc.
- print every n-th model
- for parallax models check if t_0_par is fixed and give warning, if not
- fits with 0 blending flux for some datasets
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data so that anomalies are not missed etc.
- add option to adjust Y scale to plot model fully
- add scipy to _check_imports() - requires siginificant code to be added to _check_imports() in order to find out if t_E prior is used
- in _parse_fit_constraints_prior() add a check if the priors are defined for fit parameters
- flux constraints for binary source models (note that for plotting it is now set to first dataset)
- allow Model.set_magnification_methods_parameters()
- triangle and trace plots - add option to plot burn-in as well
- methods - if only single string is provided, then this is a default method
- move _get_weighted_percentile() to a separate file with utils because it doesnt depend on self; maybe there are other similar functions
- allow LD parameters to be fitted
