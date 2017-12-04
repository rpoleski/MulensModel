# Goals for the MulensModel Package

## MulensModel Requirements

### Lensing Effects

|Effect|Implemented?|
|------|------|
| source trajectories including| |
| * annual parallax| Yes |
| * satellite parallax effects| Yes |
| point-source/point-lens magnification curves including| |
| * finite source with limb darkening| Yes, but very slow |
| binary lens magnification curves including| |
| * finite source with limb darkening| |
| * VBBL| Yes |
| * Dominik's Adapative Contouring| No |
| binary lenses caustics can be calculated and plotted together with source trajectory.| Yes |
|event properties (tE, thetaE, rE) can be calculated from physical properties of the source and lens systems (e.g. M, pi\_s, pi\_l, and mu)| Yes |
| The code contains
| * examples that show how the main high-level functions can be used| Yes |
| * tutorial | No |


### Coding Style
- The code is written in object oriented fashion
- Unit tests are presented and passed for the main high-level
   functions except the plotting ones.
- The code defines (and follows) conventions for microlensing parameters

## MulensModel Desired Properties

### Lensing Effects

|Effect|Implemented?|
|------|------|
| xallarap| No |
| automatic determination of appropriate binary lens magnification  calculation (i.e. whether or not a given approximation gives  sufficient precision)| No |
| average magnification over exposure time during caustic crossing| No |
| Can be downloaded through pip install or other standard python  package installer.| No |
| Tested on multiple operating systems (Windows, Mac, Linux)| |

### Coding Style
- Further development is possible that will include triple lens microlensing.


## Use Cases

|Num|Description|Implemented?|
|------|------|------|
|01| simple plot | Yes |
|03| model based on physical parameters | Yes|
|05| satellite parallax setup| No |
|06| WFIRST data | No |
|07| triple lenses | No |
|08| plotting | Yes |
|10| get chi2 | |
|11| clean bad data | No |
|12| simple fit PSPL | Yes |
|13| adding parallax to a model| Yes |
|14| coordinate system | No |
|15| emcee simple PSPL | Yes |
|16| raddec | Yes |
|17| magnitudes | Yes |
|18| binary equation | UC not finished! |
|19| limb\_darkening | Yes, but lacks data|
|20| binary lens instantaneous orbital motion | No |
|21| binary source modeling | No |
|22| parameters accessed as a vector | No |
|23| ModelParameters dictionary | No |
|24| Jacobian | No |

