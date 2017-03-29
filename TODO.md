#Goals for the MulensModel Package

##MulensModel Requirements

###Lensing Effects

|Effect|Implemented?|
| source trajectories including||
| * annual parallax||
| * satellite parallax effects||
| point-source/point-lens magnification curves including||
| * finite source with limb darkening||
| point-source/binary lens magnification curves including||
| * finite source with limb darkening||
| * VBBL||
| * Dominik Adpative Contouring||
| binary lenses caustics can be calculated and plotted together with source trajectory.||
|event properties (tE, thetaE, rE) can be calculated from physical properties of the source and lens systems (e.g. M, pi_s, pi_l, and mu)||
| The code contains
| * examples and tutorial that show how the main high-level functions
    can be used.||


###Coding Style
- The code is written in object oriented fashion
- Unit tests are presented and passed for the main high-level
   functions except the plotting ones.
- The code defines (and follows) conventions for microlensing parameters

##MulensModel Desired Properties

###Lensing Effects

|Effect|Implemented?|
| xallarap||
| automatic determination of appropriate binary lens magnification
  calculation (i.e. whether or not a given approximation gives
  sufficient precision)||
| average magnification over exposure time during caustic crossing||
| Can be downloaded through pip install or other standard python
  package installer.||
| Tested on multiple operating systems (Windows, Mac, Linux)||

###Coding Style
- Further development is possible that will include triple lens microlensing.


