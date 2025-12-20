MulensModel incorporates many algorithms from other work, which deserve acknowledgement in scientific publications (see [Muna et al. 2016](https://arxiv.org/abs/1610.03159)). When writing a paper that is based on MulensModel please cite MulensModel paper and also the papers for the relevant algorithms:


### Finite Source Magnification for Point lenses

`finite_source_uniform_Gould94` OR `finite_source_uniform_Gould94_direct`:
[Gould 1994 ApJ, 421L, 71](https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract)

`finite_source_uniform_WittMao94`: 
[Witt and Mao 1994 ApJ, 430, 505](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract)

`finite_source_LD_WittMao94`: 
[Witt and Mao 1994 ApJ, 430, 505](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract) and
[Bozza et al. 2018 MNRAS 479, 5157](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract)

`finite_source_LD_Yoo04` OR `finite_source_LD_Yoo04_direct`:
[Yoo et al. 2004 ApJ, 603, 139](https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract)

`finite_source_uniform_Lee09` OR `finite_source_LD_Lee09`:
[Lee et al. 2009 ApJ, 695, 200](https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract)

#### Point Lens with Shear and Convergence:

[Chang and Refsdal 1979 Nature, 282, 561](https://ui.adsabs.harvard.edu/abs/1979Natur.282..561C/abstract)

### Binary Lenses

#### Magnification Algorithms:

For all binary models (root solver): 
Skowron & Gould 2012 [ASCL](http://ascl.net/1212.005) [arXiv:1203.1034](https://ui.adsabs.harvard.edu/abs/2012arXiv1203.1034S/abstract) and 
[Bozza 2010 MNRAS, 408, 2188](https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract)

For binary models with external shear:
[Peirson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...927...24P/abstract) and [Vedantham et al. 2017 ApJ 845, 89](https://ui.adsabs.harvard.edu/abs/2017ApJ...845...89V/abstract)

quadrupole approximation: 
[Gould 2008 ApJ, 681, 1593](https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract) and 
[Pejcha & Heyrovsky 2009 ApJ, 690, 1772](https://ui.adsabs.harvard.edu/abs/2009ApJ...690.1772P/abstract)

hexadecapole approximation: 
[Gould 2008 ApJ, 681, 1593](https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract) and 
[Pejcha & Heyrovsky 2009 ApJ, 690, 1772](https://ui.adsabs.harvard.edu/abs/2009ApJ...690.1772P/abstract)

VBM/VBBL: 
[Bozza 2010 MNRAS, 408, 2188](https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract) and 
[Bozza et al. 2018 MNRAS 479, 5157](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract)
[Bozza et al. 2025 A&A 694, 219](https://ui.adsabs.harvard.edu/abs/2025A%26A...694A.219B/abstract)

Adaptive Contouring:
[Dominik 2007 MNRAS, 377, 1679](https://ui.adsabs.harvard.edu/abs/2007MNRAS.377.1679D/abstract)

#### Non-standard Parameterizations:

Cassan08:
[Cassan 2008 A&A, 491, 587](https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract)

#### Binary Lens with Shear and Convergence:

[Peirson et al. 2022 ApJ, 927, 24](https://ui.adsabs.harvard.edu/abs/2022ApJ...927...24P/abstract)

#### Keplerian orbit parameters calculation:

[Skowron et al. 2011 ApJ, 738, 87](https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract)


### Microlensing Parallax:

For microlensing parallax calculations MulensModel calls [Astropy](https://www.astropy.org/index.html) package: 
[The Astropy Collaboration et al. 2013 A&A 558, 33](https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A/abstract) and 
[The Astropy Collaboration et al. 2018 AJ 156, 123](https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract)


### Xallarap Parametrization:

[Zhai, Poleski, Zang et al. 2024 AJ 167, 162](https://ui.adsabs.harvard.edu/abs/2024AJ....167..162Z/abstract)

### Fitting Algorithms:

Our examples use EMCEE and pyMultiNest for parameter estimation. If your code uses them, then please cite the papers that presented them:

EMCEE:
[Foreman-Mackey, Hogg, Lang & Goodman 2013 PASP 125, 306](https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract)

pyMultiNest:
[Buchner et al. 2014 A&A 564, 125](https://ui.adsabs.harvard.edu/abs/2014A%26A...564A.125B/abstract).
Note that algorithms and implementation for underlying fortran code were presented by:
[Feroz & Hobson 2008 MNRAS 384, 449](https://ui.adsabs.harvard.edu/abs/2014A%26A...564A.125B/abstract), 
[Feroz, Hobson & Bridges 2009 MNRAS 298, 1601](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1601F/abstract), and 
[Feroz, Hobson, Cameron & Pettitt 2019 Open Journal of Astrophysics 2, 10](https://ui.adsabs.harvard.edu/abs/2019OJAp....2E..10F/abstract).

