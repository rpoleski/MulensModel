# MulensModel
Microlensing Modeling package

We present a python package for modeling of microlensing events. The package aims to support a number of higher order effects: microlensing parallax, finite source, binary lens (with lens rotation), and triple lenses in future.

### Second order effects

|Description|Affects<br> source<br> trajectory?|Affects<br> mag. map?|Affects<br> transformation<br> of mag. curve to<br> light curve?|   
|---|:---:|:---:|:---:|  
|flux constrained      | | | + |
|finite source         | | + | |
|limb darkening        | | + | |
|ellipsoidal source    | | + | |
|parallax              | + | | |
|xallarap              | + | | |
|multiple sources      | | | + |
|binary lens           | | + | |
|binary lens ds/dt     | | + | |
|binary lens d$\alpha$/dt | + | | |
|triple lens - future  | | + | |

**Astrometric microlensing** - we don't plan to support astrometric microlensing at this point, but we recognize that it will be used routinely in WFIRST data. Hence, we want to write the code that could be relatively easily upgraded to handle astrometric microlensing. If you see a point that could be problematic if astrometry is added, _please report it here_. 


### Model conventions

* Gould 2000 [http://adsabs.harvard.edu/abs/2000ApJ...542..785G](http://adsabs.harvard.edu/abs/2000ApJ...542..785G)
* Skowron et al. 2011 [http://adsabs.harvard.edu/abs/2011ApJ...738...87S](http://adsabs.harvard.edu/abs/2011ApJ...738...87S) - appendix A

Definitions of microlensing parameters:

* t\_eff - Gould 2013 [http://arxiv.org/abs/1312.6692](http://arxiv.org/abs/1312.6692)
* caustic crossing parameters - Cassan 2008 [http://adsabs.harvard.edu/abs/2008A%26A...491..587C](http://adsabs.harvard.edu/abs/2008A%26A...491..587C)
* single caustic crossing parameters - Albrow et al. 1999c [http://adsabs.harvard.edu/abs/1999ApJ...522.1022A](http://adsabs.harvard.edu/abs/1999ApJ...522.1022A)
* caustic size _w_ [http://adsabs.harvard.edu/abs/2009ApJ...698.1826D](http://adsabs.harvard.edu/abs/2009ApJ...698.1826D) refers to [http://adsabs.harvard.edu/abs/2005ApJ...630..535C](http://adsabs.harvard.edu/abs/2005ApJ...630..535C)
* check if new parameters are defined here: Liebig, D'Ago, Bozza, and Dominik 2015 [http://adsabs.harvard.edu/abs/2015MNRAS.450.1565L](http://adsabs.harvard.edu/abs/2015MNRAS.450.1565L)
* MORE TO BE ADDED


### How do we develop the code?

For tasks to be done see [developers board file](developers_board.md).  

How to work on new code?  
1. Write use cases if they don't exist yet.  
2. Find the objects, some of them may not be visible in use cases.  
3. Organize the objects and how they interact, if needed update [documents/Classes.md](documents/Classes.md) file.  
4. Write a few unit tests __before__ you write any function.  
5. Write the function and work on it until unit tests are passed.  

We try to follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) conventions.  

#### Naming conventions
* variables: underscore\_separated
* constants: UPPERCASE\_AND\_UNDERSCORE
* functions: underscore\_separated
  * function parameters: underscore\_separated
* class names: CamelCaseFirstLetterCapitalized
  * properties: underscore\_separated
  * private properties: \_leadingUnderscoreAndCamelCase
  * methods: underscore\_separated

### Documentation

We're using sphinx to produce documentation. To update documentation, go to `docs/` and run:
```
sphinx-build -b html source .
```

### Version numbers

Version numbers are according to MAJOR.MINOR.PATCH scheme - see [Semantic Versioning](http://semver.org/). In short:

* Patch version must be incremented if only backwards compatible bug fixes are introduced. A bug fix is defined as an internal change that fixes incorrect behavior.
* Minor version must be incremented if new, backwards compatible functionality is introduced to the public API. It must be incremented if any public API functionality is marked as deprecated. It may be incremented if substantial new functionality or improvements are introduced within the private code. It may patch level changes.
* Major version must be incremented if any backwards incompatible changes are introduced to the public API. It may include minor and patch level changes.
* Reset patch and minor version when major version is incremented. Reset patch when minor version is incremented.

