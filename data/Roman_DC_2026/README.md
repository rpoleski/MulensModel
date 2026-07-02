# Roman ephemeris - Data Challenge 2026

This directory provides information about Roman ephemeris used in the Data Challenge 2026.

## Beginner tier

The data for beginner tier contain X,Y,Z positions of Roman. The conventions are not specified, but I've figured out these are cartesian ecliptic coordinates in AU with Solar System barycenter at the origin. MulensModel allows reading the ephemeris either in specific format from JPL Horizons email format or cartesian equatorial coordinates in AU with Earth at the origin. The file `get_ephemeris.py` provides the conversion. It provides the conversion that is printed to six `RMDC26_ephemeris_season_*.dat` files. These files are also concatenated to `RMDC26_ephemeris_all_seasons.dat`. You can then read these files using `ephemerides_file` keyword of the `MulensData` class for proper inclusion of the satellite position in calculations (in models with microlensing parallax).

**Note** - the differences in magnification between Earth and Roman are at the negligible level. The situation will be different if there was another observatory taking data at the same time.

