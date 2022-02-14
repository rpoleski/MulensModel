It is important to give credit for foundations of scientific code used. When writing a paper that is based on MulensModel please cite MulensModel paper and also the papers for the relevant algorithms:


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


### Binary Lenses

#### Magnification Algorithms:

For all binary models (root solver): 
Skowron & Gould 2012 [ASCL](http://ascl.net/1212.005) [arXiv:1203.1034](https://ui.adsabs.harvard.edu/abs/2012arXiv1203.1034S/abstract) and 
[Bozza 2010 MNRAS, 408, 2188](https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract)

quadrupole and hexadecapole approximations:
[Gould 2008 ApJ, 681, 1593](https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract) and 
[Pejcha & Heyrovsky 2009 ApJ, 690, 1772](https://ui.adsabs.harvard.edu/abs/2009ApJ...690.1772P/abstract)

VBBL: 
[Bozza 2010 MNRAS, 408, 2188](https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract) and 
[Bozza et al. 2018 MNRAS 479, 5157](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract)

Adaptive Contouring:
[Dominik 2007 MNRAS, 377, 1679](https://ui.adsabs.harvard.edu/abs/2007MNRAS.377.1679D/abstract)

#### Non-standard Parameterizations:

Cassan08:
[Cassan 2008 A&A, 491, 587](https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract)


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

