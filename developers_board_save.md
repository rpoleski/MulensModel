## Specific tasks to be performed
(__boldfaced__ correspond to this month goals)

* __annual parallax calculation - test accuracy__
* __satellite parallax__
* move long unit tests to the end
* t\_0\_par is MulensTime instance
* one class per file
* one test file per class
* add check on astropy version minimum 1.2 in MulensData
* write @property for Model that returns Galactic and ecliptic coordinates based on \_coords
* MulensTime._get_date_zeropoint() and MulensData._get_date_zeropoint() - use dictionary 
* remove EMPTY files
* pass datasets from Event to Model or vice versa
* check longest files - does every function have a description?
* add a check (and warning if found) that data specified are before 1992 or after 2050
* change "type(val) is SomeType" to "isinstance(val, SomeType)"
* better import of the module so that all main classes are accesiable
* no unit tests for private functions: \_fun()
* Fit() should use marginalized distributions of fluxes
* get rid off get_jd_zeropoint from MulensData and its tests
* in unit tests if you want to assert that exception was raised then use [these](http://stackoverflow.com/questions/129507/how-do-you-test-that-a-python-function-throws-an-exception) methods
* all "#" comments are sentences
* use case 16 - code all coords features
* add Reduced JD and RHJD
* in Model: Fix ra and dec setters
* in Model: allow parameters to be set as in use case 01
* Reconsider implementation of plotting in use case 08

