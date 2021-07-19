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
- [DONE] Plotting for 1L2S, multiple datasets

## Clean Up Code
- Delete code that has been commented out:
    - [DONE] in unit tests for Event [DONE], FitData [DONE], Model [DONE]
    - [DONE]in Event[DONE], Model[DONE], MulensData[DONE], FitData[DONE]
- [DONE]Add Exceptions for methods and attributes that no longer work
- [DONE] Make usage of source_flux and f_source consistent (also for blend)
    - [DONE]Find an example to check that output of FitData can be input to
      model plotting
        - [DONE] Now, make it work properly...
    - [DONE] Add source_flux, blend_flux where necessary. Add deprecation warnings and 
      fixes for f_source, f_blend.
    - [DONE] q_flux --> flux_ratio
- [DONE] put deprecated keywords back in documentation with comment about deprecation 
and alternatives.
- [DONE] Remove Fit class (replaced with a deprecation error)

## Check All Examples Work
I think I did this, but I also might have been side-tracked by adding unit tests
- 01--16
    - [COMPLETED]: 01, 02, 03, 04, 05, 06, 07[RP], 09, 10, 11, 12, 13, 14, 15, 17
    - example_08 has not been checked (need to figure out how to open a jupyter 
      notebook)
    - example 16 is in the master branch, so can't be tested yet.
      
## Check All Previously Functioning Use Cases Work
- [COMPLETED]: 
    - working: 01, 03, 08, 10, 12, 13, 15, 16, 17, 18, 20, 23, 24, 25, 26b, 
        27, 29, 30, 31, 32, 33, 34, 37
    - not implemented: 05, 05b, 06, 07, 11, 14, 19, 22, 26, 35, 36, 38, 39
    - 21: raises NotImplemented error, but may be at least partially 
        implemented at this point.
        
- [DONE] Check with KMT about 180003 data.
    
## [DONE] Check PEP8 Compliance

# Post-merge

## Update to Version 2.0.0

## Fix alpha convention

## Check all Examples and Use Cases STILL Work
- Unit Tests completed
- Examples completed EXCEPT: 07, 08
- Example 16 appears to work
- Use cases: 01, 03, 08, 10, 12, 13, 15, 16, 17, 20, 23, 24, 25, 26b, 27, 29, 30,
    32, 33, 34, 37
  - open issues: 18
  - 31 is not implemented, contrary to above checklist.

## Check PEP8 Compliance

## Check documentation

## Write release notes explaining changes
