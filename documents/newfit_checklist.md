# Pre-merge Checklist

## Unit tests
### Event
- [DONE] tests for various combinations of fixed blending (DONE), fixed source (DONE), 
fixed q_flux (DONE)
- [DONE] tests for event.chi2 (does not change unless fit_fluxes or get_chi2 is run)


### FitData
- get_residuals() for diff values of phot_fmt
- gradient unit tests for parallax, diff combinations of u_0, t_eff, t_E (only 
do if easy)

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