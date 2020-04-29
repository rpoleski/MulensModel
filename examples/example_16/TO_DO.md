## List of task to be done:

- coordinates - to allow parallax fits
- MulensData - just provide *str*
- datasets - guessing 245/246
- data here for tests, so that PATH doesnt change
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- limit epochs in "best model" plot
- limit Y axis in "best model" plot
- trace plot
- methods for Model
- self._plots - check what is there
- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???
- should version be printed on output?
- allow plotting many models from posterior
- MulensData() - use try/except
- Add ln_prior values to blob? At some point we will want to save that information in output files
- settings['input_file_root'] = input_file_root - in final function and use it for default output files names
- check if output files (including plots) exists at the begin
- Fitting method to be added: scipy.optimize, pymultinest, ???
- add check if 't_0' is covered by data and give warning if not

