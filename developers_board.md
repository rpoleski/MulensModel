# Next big tasks - think about the order:
2. full Keplerian motion of the lens (because 2-parameter approximation is unphysical and is already implemented)
3. triple sources
4. triple lenses
5. astrometric microlensing
6. better support of FFP calculations (see "next smaller tasks" below)
7. terrestial parallax

# Next smaller tasks:
2. ModelParameters - note that t_0_1 etc. are for "binary source" models
3. make directory with reference plots from examples
4. remove unused branches
5. UC40 - finish it and implement Event.get\_lc and plot\_lc
7. UC38 - finish it
9. faster FFP calculations
11. "import MulensModel as mm" everywhere: use_cases/->30,31,34 examples/run_time_tests/check_system.py examples/checks.py TUTORIALS
12. satellite positions read from a file similar to data challenge files
13. check open issues on github
14. "add a list of public datasets:" - search below
15. In longest files re-order functions to have documentation in expected order: modelparameters.py, model.py, event.py, mulensdata.py, fitdata.py (in that order)
16. FitData.get_residuals - test for binary source

## Specific tasks to be performed
**boldfaced** tasks are most important because requested by the users

_italics_ mark important tasks

Changes for planned v3 are here: [documents/MM\_v3.md](documents/MM_v3.md)

* Install
  * **PYPI website** - some links: [1](https://docs.python.org/3.7/extending/index.html) [2](https://github.com/dariomanesku/cmft/issues/28)
  * **test Windows installation**
  * in setup.py in setup() add keywords: long\_description, classifiers
  * virtualenv; pip install -r requirements.txt; its best to install the dependencies first
  * more metadata in setup.py
  * more on setup.py: [link](https://github.com/kennethreitz/setup.py)
  * compile with "pedantic" flags for compilers
* Documentation
  * magnification\_methods.pdf - add full references there
  * Sagan workshop hands-on activity in MM
  * examples as ipython notebooks
  * Add \_\_repr\_\_ functions to Lens and Source
  * Include full documentation via setup.py data\_files mechanism.
  * note that all plotting functions require plt.show() or plt.save()
  * try removing Attributes from docstrings - just make short @property functions
  * add a note that pi\_E is "geocentric" (and "heliocentric" has the same length of vector but is rotated)
  * _example 8 corrections - PSBL, not PSPL; clarify removing the anomaly_
  * make sure that website shows correct version of MM
  * note that we are not checking for negative source or blending flux
  * add a list of public datasets: [VVV paper](https://ui.adsabs.harvard.edu/abs/2019arXiv190704339N/abstract) was published?, [second VVV paper](https://arxiv.org/abs/2106.15617) was published?, add link to K2/MCPM?; LINK the file [documents/public\_data\_list.md](documents/public_data_list.md) somewhere
* Effects:
  * **Binary source - see documents/binary\_source\_notes.md**:
    * _extract flux ratio for binary source models when fluxes are fitted via regression_
    * finish use cases
    * Fit.fit\_fluxes docstring to be updated
    * which\_parameters() - note that it doesnt work for binary source parameters, but the parameters work properly; just BSPL and rho\_2 etc. are optional
    * parallax models
    * different t\_E for each source (correct Model.set\_times)
    * test binary source with exactly one rho\_X defined
    * add t\_eff\_1, t\_eff\_2
  * **Magnification function provided by the user - already started in user\_method branch; also this could be used to model variable source events - note that**
  * **triple lens** 
    * test [ob07349](https://ui.adsabs.harvard.edu/#abs/2016AJ....152..125B/abstract) to see if center of mass can reference point for t0 and u0
    * caustics calculations
    * plan use cases and unit tests - note in [documents/TRIPLE\_LENS.md](documents/TRIPLE_LENS.md)
    * "reset" triple\_lens branch? (i.e., copy code, remove branch, start new branch e.g. tripleLens, use the code copied at begin)
    * [documents/TRIPLE\_LENS.md](documents/TRIPLE_LENS.md) - make notes in order there
    * check triple lens use cases in master branch
    * use cases, point source, hexadecapole...
  * triple source calculations
  * Finite Source
    * FSPL with low magnification - do [Witt & Mao 94](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract) or [Witt 95](https://ui.adsabs.harvard.edu/abs/1995ApJ...449...42W/abstract) give the right formulas?
    * FSPL ray shooting (ala getmag\_rs\_single.f)
    * Yoo+04 full formalism 
    * get gamma/u LD coeffs from Claret papers etc.
    * [Lee+09](https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract) - gradient calculations for uniform source, also faster calculations - profile
    * FSPL for large sources using [Agol 2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...594..449A/abstract)
  * Quadratic limb darkening
  * Multi-lens ray shooting:
    * mapmaking version which adds new rays as needed (but remember that it runs for fixed (s,q) only!)
    * Yossi idea to find all the images
  * Orbital motion like in [VBBL 2.0](https://arxiv.org/abs/1805.05653)
  * calculate jerk parallax degeneracy: [Park+04](https://ui.adsabs.harvard.edu/abs/2004ApJ...609..166P/abstract) [Gould 04](https://ui.adsabs.harvard.edu/abs/2004ApJ...606..319G/abstract)  
  * topocentric/Earth parallax
  * Chang-Refsdal binary calculations
  * elliptical source magnification [Heyrovsky & Loeb 1997](https://ui.adsabs.harvard.edu/abs/1997ApJ...490...38H/abstract)
  * magnification calculated for a set of points, not just a trajectory - this way we could, e.g., plot magnification maps
* Parameterization
  * Cassan 2008 binary lens parameters:
    * option to change scaling (from [0,1] to C08 params) to work well near topology change
  * [Albrow et al. 1999](https://ui.adsabs.harvard.edu/abs/1999ApJ...522.1022A/abstract) (also Cassan 2008 Sec. 5)
  * t\_eff as a parameter for large u\_0 - see [Andy paper](https://arxiv.org/abs/1312.6692) and maybe also other from [Jen 2012 paper](https://ui.adsabs.harvard.edu/abs/2012ApJ...755..102Y/abstract), i.e., f\_lim=f\_s/u\_0 and q\*t\_E
  * caustic size w [Dong+09](https://ui.adsabs.harvard.edu/abs/2009ApJ...698.1826D/abstract) refers to [Chung+05](https://ui.adsabs.harvard.edu/abs/2005ApJ...630..535C/abstract); see also Skowron+11 Eq.2 
  * check if new parameters are defined here: [Liebig, DAgo, Bozza, and Dominik 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450.1565L/abstract)
  * [Heyrovsky 2003](https://ui.adsabs.harvard.edu/abs/2003ApJ...594..464H/abstract) parametrization of limb-darkening
  * t\_0\_planet, u\_0\_planet, t\_E\_planet instead of s, q, alpha
  * [Dominik 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.393..816D/abstract) for PSPL
  * [Jung+17](https://ui.adsabs.harvard.edu/abs/2017AJ....153..129J/abstract) - rotating triple lens - somehow special version of xallarap
* Function Improvements/Expansion:
  * BinaryLens class:
    * _VBBL2.0 - are we using accuracy limit as default? If so then we should switch to relative accuracy_
    * should BinaryLens() accept source\_x/y as lists or arrays?
    * function for center of mass shift (currently: shift\_x in trajectory.py, x\_shift in binarylens.py, xcm\_offset in caustics.py)
    * central and planetary caustic properties: [Chung et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...630..535C/abstract) and [Han 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...638.1080H/abstract)
    * consider using Utils.complex\_fsum() in BinaryLens functions: \_jacobian\_determinant\_ok\_WM95()
    * faster hexadecapole using [Cassan 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.468.3993C/abstract) ([code](https://github.com/ArnaudCassan/microlensing/blob/master/microlensing/multipoles.py))
    * there should be an option for the user to ignore "wrong number of solutions" error and replace it with a warning or dismiss it fully
    * error message for not 3 or 5 images - This can be improved by checking the denominators in the equation (before polynomial is fully formed - see Zoey Mathematica notes)
    * error message for not 3 or 5 images - in "Consider ..." note that one must use Model.set\_magnification\_methods() for epochs far from peak
    * VBBL additional parameter so that it works for huge source and tiny caustic
  * Caustics class:
    * Caustics.\_calculate - optimize using vectors instead of a loop
    * _Caustics calculations using [Erdl & Schneider 1993](https://ui.adsabs.harvard.edu/abs/1993A%26A...268..453E/abstract) approach; there is also [Asada 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...390L..11A/abstract) and I dont know if it makes sense to code that one as well_
    * root solver can be used the same as in binarylens.py - not needed for binary lens and Erdl & Schneider 1993
    * smaller points
    * correct conditions in get\_caustics()
  * Event class:
    * **Allow fluxes to be fixed in chi^2 calculation (e.g. given a particular fs, fb, which you might want to do if you want fs as a chain parameter); also think how it will work for binary sources**
    * **give access to all fluxes without changing data\_ref**
    * **plot magnitude difference between 2 models for residuals plot**
    * _it seems it doesnt have plot\_trajectory()_
    * add plot\_source() etc. already implemented in Model
    * Event should sync information on which of the 3 types of parallax are used, so that if it is specified for event, then there will be exception if one dataset is missing earth\_coords etc. In general there should be some way to make sure which parallax types are used in which calculation of magnification.
    * Class Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
    * reduce calls to Fit.fit\_fluxes()
    * chi2\_gradient() should cope NaN values in a way similar to get\_chi2()
    * **check all functions that should pass** fit\_blending parameter - Event.plot\_model, what else??? Already done: Event.get\_ref\_fluxes()
    * chi2 with maximum value provided - if the chi2 for point-source gives chi2 larger than specified limit, then finite source calculations are not undertaken (this should significantly speed-up MultiNest)
    * get flux and its error in reference system
    * change order to improve the website
    * gradient - fluxes as well? if so, then start using the second test in test\_event\_chi2\_gradient()
    * for consistency, it would be good to combine get\_chi2\_for\_dataset() and get\_chi2\_per\_point()
    * other likelihoods, e.g., [Dominik+18](https://arxiv.org/abs/1808.03149), [ARTEMiS](https://ui.adsabs.harvard.edu/abs/2008AN....329..248D/abstract), SIGNALMEN, [RoboTAP](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A..55H/abstract)
    * function that calculates cumulative chi2 so that it can be plotted easier
    * Binary source - see optimization comment at begin of Event.get\_chi2\_for\_dataset()
    * plot cumulative chi2 difference between 2 models - magnification and magnitude spaces
  * Fit class:
    * should use marginalized distributions of fluxes (if those are from linear fits); JCY - it needs UC
    * n\_sources somehow inconsistent in different places
  * Horizons class:
    * JPL Horizons
      * correct JPL Horizons => CSV file format; also example usage
      * check if Horizons e-mail is for correct satellite
    * BJD
      * conversions to BJD from HJD, JD etc. ([astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections))
      * BJD\_TDB in satellite ephemeris [astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections)
  * Lens class:
    * \_\_repr\_\_ function needs work                                         
    * a\_proj, couples with source distance in mulensmodel to determine s.  
    * 2-body example 3 is missing s. Why? Does that work?                  
    * problem with tracking number of masses, esp when successively defining masses (see test\_Lens.py)
    * implement triple+ systems  
  * MagnificationCurve class:
    * what to do if some magnifications are nan, inf, or None?; maybe give warning and return None in these places
    * re-write magnification() to use lazy loading (here or in model.py)
  * Model class:
    * **in functions magnification(), plot\_magnification(), and plot\_trajectory() use satellite\_skycoord from \_\_init\_\_ if available**
    * **plot\_lc() - add satellite option like in plot\_magnification(), other options as well - use keywords passed to self.magnification()**
    * reorder functions so that it looks good on website
    * bandpass option could simultaneously account for limb-darkening and source_flux_ratio for a given band (to be applied in a few options that have bandpass option)
    * Model.set\_parameters() should remember previously set values (of course unless they are overwritten)
    * Class Model should not allow accessing attributes that should not be there, eg., q for single lens case.
    * Function that prints RA, Dec, t\_0\_par, t\_0\_kep, types of parallaxes turned on, satellite info, limb coeffs and textual description of type of model
    * plot\_trajectory() - mark epochs using colorscale? Maybe it should be passed by kwargs (if so, then add example)
    * Should get\_satellite\_coords() use caching?
    * we should have versions of all plot functions to use magnifications instead of magnitudes; also add access via Event
    * add option to use random zorder when plotting multiple datasets (e.g. gaussian with sigma depending on number of plotted datapoints)
    * does Model.plot\_lc() give docstring for flux\_ratio\_constraint option? If not, then present it in plot\_magnification() docstring
    * use plt.axhline() to plot 0 line in residuals plots at the end, using t\_min,t\_max
    * get\_residuals needs unit test
    * plot\_trajectory - mark data epochs (as pyplot kwargs use MulensData.plot\_properties)
    * plot\_data & plot\_residuals - if color is not set by the user and show\_bad is True, then X-s are plotted using different color
    * allow rotating plots, so that source travel "exactly" from left to right
    * when plotting residuals allow keywords to be passed for plotting 0-line
    * get\_data\_magnification() - should there be data.good masking?
    * get\_ref\_fluxes - use caching (self.\_fit) inside it
    * plot\_lc() should have satellite\_skycoord keyword as plot\_magnification() has
    * plot\_source and binary sources - allow different kwargs for each of the sources?
    * check\_if\_caustic\_crossing() function can be added - for non-parallax, no-orbital motion it requires calculation of cusps positions for each caustic and checking if they are on the same side of the trajectory
    * try/except in pyplot commands and catch AttributeError - these may be miss-spelling etc.
  * ModelParameters class:
    * _values in dimensionless astropy.quantity should be changed to float, other types should be rejected (unless it is a time unit etc.)_
    * _LaTeX strings with parameters names (useful e.g. for corner plots or tables with results)_
    * a check if each PARAMETER.setter function changes parameter that was previously set
    * Transform t\_E and other parameters between geocentric and heliocentric frames.
    * option to return alpha, dalpha\_dt, and ds\_dt as floats instead of astropy.quantities
    * why .rho returns None if it is not defined? In other similar cases we have KeyError. Should that be changed? (if so, then maybe only after changing version to 2.0.0)
    * to make \_check\_valid\_combination\_1\_source\_...() shorter, make a boolean dict that says if given parameter is defined or not
    * change order to improve the website
    * check if t\_eff and t\_star can be used as input simultaneously
    * check if input values are floats (or something else accepted)
    * warning if t\_0\_par is changed too many times
    * add is\_Cassan08 or is\_standard\_binary?
    * \_\_repr\_\_ - if Cassan 2008 then should we print standard parameters as well?
  * MulensData class:
    * **Errorbar scaling, in particular the two parameter.**
    * _quick look alignment for MulensData objects - just use interpolation_
    * make \_\_init\_\_() shorter
    * add version of n\_epochs that uses only good epochs
    * read settings from file header: flux vs. mag, filter, satellite info
    * change order to improve the website
    * docstring phot\_fmt vs. input\_fmt
    * data\_and\_err\_in\_input\_fmt() and Fit.get\_input\_format() possible can be deprecated or removed because we shifted to chi2\_fmt instead of input\_fmt
    * when plotting data, make sure that max/min limits on Y axis include errorbars, if the errorbars are shown
    * export/save given data file in scale of other dataset and model
    * this line may be wrong for some values of char: kwargs['fmt'] = kwargs['fmt'].replace(char, "")
    * for plt.scatter() the color can be set as 'facecolor', 'facecolors', or 'edgecolors' and this should be dealt with in \_set\_plot\_properties()
    * for plotting X for bad data use large size and/or thinner line
    * separate colors (or even kwargs) for X-es as an option (to get contrasting colors see https://itsphbytes.wordpress.com/2016/08/29/complementary-colors-python-code/)
  * PointLens class:
    * in \_integrand\_Lee09\_v2() use loop and math.fsum in "values = ...", check performance as u\_ and theta\_ have shape of (90, 1000)
    * make WM method faster: 1) interpolation done once for many values; 2) interpolate different function; 3) allow changing number of annuli; 4) divide A by 2 different functions (z<1 or z>1) and interpolate these - like in VBBL (arguments of interpolation are: 1) (rho, z), 2) (log rho, log z) etc.)
    * add [Witt and Atrio-Barandela 2019](https://arxiv.org/abs/1906.08378)?
    * get\_pspl\_magnification() - change it to operate on u^2, not u, so that np.sqrt() calls are reduced
    * 1+2/u^4 approximation for very large u
    * try to remove sympy and use scipy instead (not possible in Sep 2020, because there is no elliptical integral of third kind in scipy)
  * SatelliteSkyCoord class:
    * attach magnification\_methods to SatelliteSkyCoord so that they overwrite Model and MagnificationCurve settings when given SatelliteSkyCoord is used
  * Trajectory class:
    * _warning when too many annual parallax calculations are conducted_
    * _\_get\_delta\_satellite() should be using self.times_
    * annual parallax caching - if moved to MulensData, then would be faster because hashing of coords and time vector takes time
    * colorscale time or magnification (see Fig 2 in [Ranc+19](https://arxiv.org/abs/1810.00014))
    * plot in units of theta\_star (or even days) instead of theta\_E
  * UniformCausticSampling class:
    * add function that calculated t\_in & t\_out in analogy to get\_x\_in\_x\_out() - the problem is that for given (u\_0, alpha) - there are 2, 4, or 6 caustic crossing points
    * prior probability (in get\_uniform\_sampling()) of hitting given caustic based on (d, q, alpha) - see [Kains+12](https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.2228K/abstract)
    * prior only in given box (in get\_uniform\_sampling())
    * prior includes t\_E prior provided by user (also add option to use Mroz+19 etc) in get\_uniform\_sampling()
    * trajectories that cross 2 different caustics
    * profile calculations
  * Utils class:
    * in np.any() ifs give more information in warning e.g., "out of 1234 values provided, the fails are: 12, 345, 678 (0-based)"
    * add u(a) function: u = np.sqrt(2A/np.sqrt(A^2-1.) - 2.)
    * utils.py:57 - code produces 2 warnings, should produce just one; use masking
    * documentation - use ITALICS
    * uses masks there if warnings and apply None etc. - check examples after change!
    * add a function that calculates mag/flux in a different flux space - we currently do it in MulensData.plot() and Model.get\_residuals()
  * MulensObjects submodule:
  * Plotting:
    * _plt.plot() .scatter() and .errorbar() should share default colors, so that when you plot one dataset and one model and dont set colors, then they are of different colors_
    * for plotting functions option to pass pyplot.Axis and pyplot.Figure instances and call e.g. Axis.scatter() instead of pyplot.scatter(); for a simple example see [here](https://github.com/rpoleski/K2-CPM/blob/master/source/K2CPM/plot_utils.py)
    * subplots with shared X-axis (plt.subplots(2, 1, sharex=True, gridspec\_kw={'height\_ratios': [4, 1]}, figsize=???, dpi=100)) - start in Example 5
    * add option to plot satellite coordinates as in Henderson+16 where K2 and Spitzer orbits were compared; i.e., make Trajectory.\_get\_delta\_satellite() public and add appropriate example
    * add plotting with fit\_blending=False for functions that use magnitude space
    * add plt.xlim() and ylim in plotting functions (using t\_start) etc.; then also update (simplify) tutorials, examples etc.
    * caustics for trajectory plot with single lens models
    * plot upper limits instead of photometry with negative flux
  * Examples and use cases:
    * use case 34 - make the output meaningful 
    * example 4 - use "with open() as ...:"
    * _Hamiltonian MCMC [link 1](http://arogozhnikov.github.io/2016/12/19/markov_chain_monte_carlo.html) and [link 2](https://theclevermachine.wordpress.com/2012/11/18/mcmc-hamiltonian-monte-carlo-a-k-a-hybrid-monte-carlo/) and [link 3](https://colindcarroll.com/2019/04/11/hamiltonian-monte-carlo-from-scratch/)_
    * example with [parallel EMCEE](https://emcee.readthedocs.io/en/stable/tutorials/parallel/)
    * _plot many models from posterior_
    * **scipy.curve\_fit() and print parameter uncertainties**
    * **corner plots; they require [corner](https://github.com/dfm/corner.py), [pyhdust](https://github.com/danmoser/pyhdust), or [pyGTC](https://pypi.org/project/pyGTC/)**
    * _F\_s for MOA data for MB08310 differs from Janczak paper - is it caused by FSPL vs. FSBL models?_
    * add example that shows 'log\_' in the name of the parameter; central caustic anomaly planet would be best,
    * _example 08: PSBL, not PSPL_
    * gaussian process - see [Li+19](https://ui.adsabs.harvard.edu/abs/2019arXiv190407718L/abstract) (maybe it requires input in Fit class, hence would not be just an example)
    * Dan FM approach for periodic variables - two variables instead of one (identified by some postfix) - make it for alpha or x\_caustic\_in/out
    * add illustration on how to remove airmass trends
    * add example of fitting PSPL model using [Albrow (2004)](https://ui.adsabs.harvard.edu/abs/2004ApJ...607..821A/abstract) method [link](https://github.com/MichaelDAlbrow/SingleLensFitter/blob/master/SingleLensFitter.py)
    * plotting - sharex where possible
    * note in PSPL tutorial about plotting data in MulensData
    * add example that after calling Event.get\_chi2() use Event.fit to get e.g. magnifications so that the magnification is not calculated twice
    * **satellite data fitted and plotted - what is missing now?**
    * some cfg files use "../data/..." - change it to MM.DATA\_PATH somehow
    * in emcee we should check if all the starting points are in prior
    * check if MM correctly interacts with scipy.optimize.leastsq and maybe add an example
    * add Coordinates.velocity\_of\_Earth() example
    * Example 13 - make x\_caustic\_in/\_out periodic variables
    * add example with well-known code like [https://mc-stan.org/](https://mc-stan.org/)
    * MulensData.bad needs example note that one has to substitute full vector, not single values
  * Miscellaneous:
    * _COVERAGE : "coverage run --source MulensModel -m py.test" and then "coverage report" or "coverage report -m" or "coverage html" (and then open htmlcov/index.html); https://coverage.readthedocs.io/en/v4.5.x/_
    * u\_0 sign for satellite or just parallax model - some way of following u(t) evolution
    * warnings - give type for each one of them
    * Dave idea on caustic-crossing epochs as parameters t\_cc1 and t\_cc2 - see david.p.bennettATnasa.gov e-mail on Feb 11, 2019
    * when checking units use Unit.physical\_type - search for physical\_type in mulensobjects/lens.py as an example; to find places to be changed search for "isinstance" (to find these places run grep isinstance \*py mulensobjects/\*py | grep Quantity
    * use lazy loading in MagnificationCurve.magnification and/or Model.magnification
    * guessing parameters of PSPL model ([Kim+17](https://arxiv.org/abs/1703.06883) as an example)
    * add calculation of Caustic Region of Influence (CROIN) - [Penny 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...790..142Y/abstract)
    * anything from use cases that does not work yet -- see TODO.md file
    * comments at begin of each use case and example
    * use case 25 doesnt work
    * OB161195 is used only in UC25 - change it and remove directory
    * interaction with fitting routines - see [list of them](https://arxiv.org/abs/1711.03329)
    * caching of results in trajectory.py should stop at some point - if the user changes t\_0\_par or coords, then there is no point in remembering huge indexes (whole self.times)
    * profile the code (python -m cProfile script.py - you only should check 2nd and 4th column)
    * Leap seconds library - [barycorrpy](https://arxiv.org/abs/1801.01634)
    * [documents/TODO.md](documents/TODO.md) file - move content and remove
    * add transformation Jacobians, see: Batista+11 and Skowron+11
    * use cython, numba, pypy or similar (numba seems best) to speed-up calculations in Lee+09 and/or Cassan08
    * remove "import *" from \_\_init\_\_.py and mulensobjects/\_\_init\_\_.py
* Other Tests:
  * test\_event\_chi2\_gradient() - add parallax without flux gradient
  * add unit tests for Horizons and MulensData.satellite\_skycoord
  * Coordinates - write tests, possibly remove test\_Coords.py
  * t\_eff is not tested
  * plt.scatter -> plt.plot; after that we can start unit tests for plt calls
  * check for memory leaks by running long calculations and monitoring RAM usage
* Style/Architecture:
  * Are we consistent with PEP8? [check here](http://pep8online.com/) or pycodestyle command
  * PEP8 for tests/ (test\_Event is already done: test\_Event test\_ModelParameters test\_Model\_Parallax
  * for examples/
  * better import of the module so that all main classes are accessible (use \_\_all\_\_ = [...] in all files?)
  * Utils - Make subpackage/submodules that group related functions (e.g. flux2mag conversions)?

### reStructuredText:
[1st tutorial](http://gisellezeno.com/tutorials/sphinx-for-python-documentation.html), 
[2nd tutorial](http://www.sphinx-doc.org/en/stable/rest.html), 
[example](https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html)

### Xallarap references:

* [Griest & Hu 1992](https://ui.adsabs.harvard.edu/abs/1992ApJ...397..362G/abstract), 
* [Han & Gould 1997](https://ui.adsabs.harvard.edu/abs/1997ApJ...480..196H/abstract), 
* [Dominik 1998](https://ui.adsabs.harvard.edu/abs/1998A%26A...329..361D/abstract), 
* ob9919 - [Smith et al. 2002](https://ui.adsabs.harvard.edu/abs/2002MNRAS.336..670S/abstract), 
* [Ghosh et al. 2004](https://ui.adsabs.harvard.edu/abs/2004ApJ...615..450G/abstract), 
* [Jiang et al. 2004](https://ui.adsabs.harvard.edu/abs/2004ApJ...617.1307J/abstract), 
* [Poindexter et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...633..914P/abstract) - 23% of events are affected by xallarap,
* [2008 conference posting](https://ui.adsabs.harvard.edu/abs/2008mmc..confE..37R/abstract),
* [Dong et al. 2009](https://ui.adsabs.harvard.edu/abs/2009ApJ...695..970D/abstract), 
* ob07368 - [Sumi et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJ...710.1641S/abstract), 
* ob07514 - [Miyake+12](https://ui.adsabs.harvard.edu/abs/2012ApJ...752...82M/abstract), 
* mb10328 - [Furusawa et al. 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...779...91F/abstract), 
* ob130911 - [Miyazaki et al. 2019](https://arxiv.org/abs/1912.09613), 
* Roman predictions - [Miyazaki+20](https://arxiv.org/abs/2010.10315),
* ob150845 = mb15277 - Calen leads, 
