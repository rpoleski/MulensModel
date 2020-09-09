# Pre-merge Checklist

## Unit tests
### Event
- [DONE] tests for various combinations of fixed blending (DONE), fixed source 
(DONE), 
fixed q_flux (DONE)
- [DONE] tests for event.chi2 (does not change unless fit_fluxes or get_chi2 is 
run)


### FitData
- [DONE] get_residuals() for diff values of phot_fmt
- [SKIP] gradient unit tests for parallax, diff combinations of u_0, t_eff, t_E (only 
do if easy)

## Implementation
### Event
- [DONE] plot_source_for_datasets()
    1) [DONE] identify appropriate example for testing --> Modify 
    example_05_MB08310 to include trajectory plot
    2) [DONE] implement

## Clean Up Code
- Delete code that has been commented out:
    - [DONE] in unit tests for Event [DONE], FitData [DONE], Model [DONE]
    - [DONE]in Event[DONE], Model[DONE], MulensData[DONE], FitData[DONE]
- [DONE]Add Exceptions for methods and attributes that no longer work
- Make usage of source_flux and f_source consistent (also for blend)
    - [DONE]Find an example to check that output of FitData can be input to
      model plotting
        - Now, make it work properly...
    - Add source_flux, blend_flux where necessary. Add deprecation warnings and 
      fixes for f_source, f_blend.
    - q_flux --> flux_ratio
- put deprecated keywords back in documentation with comment about deprecation 
and alternatives.
- Remove Fit class

## Check All Examples Work
I think I did this, but I also might have been side-tracked by adding unit tests
- 01--16
    - Completed: None

## Check PEP8 Compliance