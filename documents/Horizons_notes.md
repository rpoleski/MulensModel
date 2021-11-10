# Notes on output format from JPL Horizons

### HORIZONS OUTPUT FILE:
  Prior to 1962, times are UT1. Dates thereafter are UTC. Any 'b' symbol in
the 1st-column denotes a B.C. date. First-column blank (" ") denotes an A.D.
date. Calendar dates prior to 1582-Oct-15 are in the Julian calendar system.
Later calendar dates are in the Gregorian system.

  Time tags refer to the same instant throughout the solar system, regardless
of where the observer is located. For example, if an observation from the
surface of another body has an output time-tag of 12:31:00 UTC, an Earth-based
time-scale, it refers to the instant on that body simultaneous to 12:31:00 UTC
on Earth.

  The Barycentric Dynamical Time scale (TDB) is used internally as defined by
the planetary equations of motion. Conversion between TDB and the selected
non-uniform UT output time-scale has not been determined for UTC times after
the next July or January 1st. The last known leap-second is used as a constant
over future intervals.

### Eastman et al. 2010 - source file gives (A1.1):
\subsubsection{HORIZONS}
JPL's HORIZONS ephemeris calculator, which is used by many to
calculate the BJDs from space telescopes, and can be used to calculate
ground-based BJDs, returns the time in \jdct = \teph = \jdtdb when the
ephemeris type is ``Vector Table''. Any conversion that uses a
HORIZONS ephemeris in \jdtdb \ but indexes it with \jdutc, as had been
done by several people we spoke with, will calculate \bjdutcp, which
can be offset from the true \bjdutc \ by up to 13.4 ms (as shown in
Fig. \ref{fig:utcvutc}), and offset from the uniform \bjdtdb \ by more
than 1 minute (as shown in Fig. \ref{fig:tdbvutc}).

### http://ssd.jpl.nasa.gov/?horizons_doc#time

 The program's time-span prompts indicate the earliest & latest dates
 that may be used for the selected target/center combination, as well
 as the type of time assumed being input (UT, TDB, or TT).

For cartesian coordinates or osculating elements tables, only TDB may
be used. For "observer tables", output may be either UT or TT. TO
CHANGE THE UT DEFAULT for observer tables, append a "TT" when entering
START time. To switch back, append a "UT" to the start time.

The three time systems are described as follows:

TDB
    ("Barycentric Dynamical Time"); typically for cartesian,
    osculating element, and close-approach tables. The uniform time
    scale and independent variable of the planetary ephemeris
    dynamical equations of motion.

TT
    ("Terrestrial (Dynamic) Time"), called TDT prior to 1991, used for
    observer quantity tables. This is proper time as measured by an
    Earth-bound observer and is directly related to atomic time, TAI.
    TT periodically differs from TDB by, at most, 0.002 seconds.

UT
    is Universal Time This can mean one of two non-uniform time-scales
    based on the rotation of the Earth. For this program, prior to
    1962, UT means UT1. After 1962, UT means UTC or "Coordinated
    Universal Time".  Future UTC leap-seconds are not known yet, so
    the closest known leap-second correction is used over future
    time-spans.

