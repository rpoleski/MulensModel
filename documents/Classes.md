## Classes description


### Class Event
* model - instance of Model
* datasets - list of UlensData instances
* fit - instance of Fit
* get_chi2()
* chi2

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
* A - array
* flux - array
* mag - array
* theta_E
* bandpass
* source_flux
* blend_flux

### Class Fit
* t - list of arrays
* mag_data - list of arrays
* mag_model - list of arrays
* res - list of arrays

### Class Lens
* set_componenet()
* n_components
* component[]
* s
* q
* m_1
* m_2
* distance
* and for future triple lens:
  * s_2
  * q_2
  * s_3
  * q_3
  * s_23
  * beta_23

### Class Source
* distance

### Class UlensData
* satellite
* t
* hjd
* mag
* err_mag


