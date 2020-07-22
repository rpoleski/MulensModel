# Pre-merge Checklist

## Unit tests
### Event
- tests for various combinations of fixed blending, fixed source, fixed q_flux
- tests for event.chi2 (does not change unless fit_fluxes or get_chi2 is run)


### FitData
- get_residuals() for diff values of phot_fmt
- gradient unit tests for parallax, diff combinations of u_0, t_eff, t_E

## Model
- test get_trajectory:
  - straight-up trajectory
  - case with annual parallax (check test_Model_Parallax.py)
  - case with satellite parallax (check test_Model_Parallax.py)
  - coords is propagating correctly (check test_Model_Parallax.py)
- test set_times:
  - keywords to test:
    t_range=None, t_start=None, t_stop=None, dt=None, n_epochs=1000
- test set_default_magnification_method:
  - change from default value
- test get_satellite_coords: (check test_Model_Parallax.py)
  - returns None if no ephemerides file set

## Implementation
### Event
- plot_source_for_datasets()
    1) identify appropriate example for testing

## Clean Up Code
- Delete code that has been commented out:
    - in unit tests for Event, FitData, Model
    - in Event, Model
- Add Exceptions for methods and attributes that no longer work
- Remove Fit class

## Check All Examples Work
I think I did this, but I also might have been side-tracked by adding unit tests
- 01--15