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

* Gould 2000 [https://ui.adsabs.harvard.edu/abs/2000ApJ...542..785G/abstract](https://ui.adsabs.harvard.edu/abs/2000ApJ...542..785G/abstract)
* Skowron et al. 2011 [https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract](https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract) - appendix A


### How do we develop the code?

For tasks to be done see [developers board file](developers_board.md).  

How to work on new code?  
1. Write use cases if they don't exist yet.  
2. Find the objects, some of them may not be visible in use cases.  
3. Organize the objects and how they interact.  
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

If the documentation is not updating, try running sphinx-apidoc again (in `docs/)`:
```
sphinx-apidoc -e -f -o source/ ../source/MulensModel/
```

To change the content of index.html, modify the file `source/index.rst`

### Version numbers

Version numbers are according to MAJOR.MINOR.PATCH scheme - see [Semantic Versioning](http://semver.org/). In short:

* Patch version must be incremented if only backwards compatible bug fixes are introduced. A bug fix is defined as an internal change that fixes incorrect behavior.
* Minor version must be incremented if new, backwards compatible functionality is introduced to the public API. It must be incremented if any public API functionality is marked as deprecated. It may be incremented if substantial new functionality or improvements are introduced within the private code. It may patch level changes.
* Major version must be incremented if any backwards incompatible changes are introduced to the public API. It may include minor and patch level changes.
* Reset patch and minor version when major version is incremented. Reset patch when minor version is incremented.

List of files to be updated: README.md docs/source/conf.py docs/source/install.rst source/MulensModel/version.py

Also for code releases you should temporarily remove use cases and developers\* files.

