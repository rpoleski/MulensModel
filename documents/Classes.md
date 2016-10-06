## Classes description


### Class Event
* model - instance of Model
* datasets - list of MulensData instances
* fit - instance of Fit
* chi2
* get_chi2()
* clean_data()
* estimate_model_params()

### Class Model
* lens - instance of Lens
* source - instance of Source
* caustic - has properties x & y
* trajectory - has properties x & y
* ra
* dec
* t_0
* u_0
* t_E
* rho
* alpha
* s - property of subclass Lens ???
* q - property of subclass Lens ???
* time - array
* hjd - array
* mag - array
* magnification - array
* flux - array
* mag - array
* theta_E
* bandpass
* source_flux
* blend_flux
* pi_E
* t_0_par
* time_data
* mag_data
* parameters()
* parallax()
* plot_lightcurve()
* plot_trajectory()
* finite_source()

### Class Fit
* time - list of arrays
* mag_data - list of arrays
* mag_model - list of arrays
* res - list of arrays

### Class Lens
* n_components
* component[]
* s
* q
* mass_1
* mass_2
* distance
* set_componenet()
* and for future triple lens:
  * s_2
  * q_2
  * s_3
  * q_3
  * s_23
  * beta_23

### Class Source
* distance
* angular_size
* I_mag

### Class MulensData
* satellite
* time
* hjd
* mag
* err_mag
* bad
* flux
* err_flux
* mag_fmt

