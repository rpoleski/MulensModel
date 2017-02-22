import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from MulensModel.utils import Utils

#There is nothing that specifies Time reference frame and scale
#(i.e. HJD_UTC vs. BJD_TDB.
# Radek: See copy-pasted texts at the bottom of this file.


class Horizons(object):
    def __init__(self, file_name=None):
        self._time = None
        self._xyz = None
        self._names = {} # This dictionary records the HORIZONS column names 
                         # (values) and normal ways we understand them (keys).
        self._names['date'] = 'Date__(UT)__HR:MN'
        self._names['ra_dec'] = 'R.A._(ICRF/J2000.0)_DEC'
        self._names['distance'] = 'delta'
        apply_float = ['distance'] # List of columns to which float() 
                                   # operator is applied.
        apply_date_change = ['date'] # List of columns to which date
                                     # format change function is applied. 
        self.data_lists = {} # Dictionary with all the data read from 
                             # HORIZONS table. Keys are the same as in 
                             # self._names and values are lists of data read. 
        for key in self._names:
            self.data_lists[key] = []
            
        self._read_horizons_file(file_name)
        
        for key in apply_float:
            for j in range(len(self.data_lists[key])):
                self.data_lists[key][j] = float(self.data_lists[key][j])
        for key in apply_date_change:
            for j in range(len(self.data_lists[key])):
                self.data_lists[key][j] = Utils.date_change(self.data_lists[key][j])

    def _read_horizons_file(self, file_name):
        """reads standard output from JPL Horizons"""
        mode_save = "START"
        with open(file_name) as in_file:
            mode = None # Other possible values are: 'header', 'header_done', 'main', and 'finished'.
            # This gives information in which part of the Horizons file we are currently in.
            for l in in_file.readlines():
                line = l[:-1]
                if len(line) == 0:
                    continue
                mode_save = mode # Value of mode variable is remembered before the line is parsed. 
                
                if line[0] == '$':
                    mode = self._process_dollar_line(line, mode)
                    if mode == 'finished': # Mode 'finished' means we have arrived at the line 
                                           # that indicates end of data table
                        break
                elif line[0] == '*':
                    mode = self._process_star_line(line, mode)
                        
                if mode != mode_save:
                    continue # Mode has changed so current line is done. 
                    
                if mode == 'header':
                    (char_beg, char_end) = self._parse_header_line(line)
                elif mode == 'main': # This is where we parse main data table.
                    for key, value in char_beg.items():
                        text = line[value:char_end[key]]
                        self.data_lists[key].append(text)
                elif mode != None:
                    raise ValueError('unexpected input: {:}'.format(line))
    
    def _process_dollar_line(self, line, mode):
        """deals with the line that starts with '$'"""
        if line == '$$SOE': # $$SOE marks start of data table.
            if mode == 'header_done':
                mode = 'main'
            else:
                raise ValueError('error: {:}'.format(line))
        elif line == '$$EOE': # $$EOE marks end of data table.
            if mode == 'main':
                mode = 'finished'
            else:
                raise ValueError('error: {:}'.format(line))
        else:
            raise ValueError('unexpected input: {:}'.format(line))        
        return mode

    def _process_star_line(self, line, mode):
        """deals with the line that starts with '*'"""
        if len(line) != 79:
            if mode is None:
                mode = 'header'
            elif mode == 'header':
                mode = 'header_done'
        else:
            if mode == 'main':
                raise ValueError('error: {:}'.format(line))
        return mode

    def _parse_header_line(self, line):
        """parse line with table header from JPL Horizons"""
        char_beg = {}
        char_end = {}
        for (key, value) in self._names.items():
            find = line.find(value)
            if find == -1:
                raise ValueError('Header key not found: {:}'.format(value))
            char_beg[key] = Utils.last_non_space_char_before(line, find) + 1
            char_end[key] = find + len(value)
        return (char_beg, char_end)

    @ property
    def time(self):
        """return time vector in whatever format was input"""
        if self._time is None:
            times = Time(self.data_lists['date'], format='iso', scale='utc') # Currently we assume HORIZONS works in UTC.
            self._time = times.tdb.jd 
        return self._time

    @ property 
    def xyz(self):
        """return X,Y,Z positions"""
        if self._xyz is None:
            self._xyz = SkyCoord(self.data_lists['ra_dec'], 
                                distance=self.data_lists['distance'], 
                                unit=(u.hourangle, u.deg, u.au)).cartesian
        return self._xyz



"""
HORIZONS OUTPUT FILE:
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
"""

"""
Eastman et al. 2010 - source file gives (A1.1):
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
"""

"""
http://ssd.jpl.nasa.gov/?horizons_doc#time
 The program's time-span prompts indicate the earliest & latest dates 
 that may be used for the selected target/center combination, as well 
 as the type of time assumed being input (UT, TDB, or TT).

For cartesian coordinates or osculating elements tables, only TDB may 
be used. For "observer tables", output may be either UT or TT. TO 
CHANGE THE UT DEFAULT for observer tables, append a "TT" when 
entering START time. To switch back, append a "UT" to the start time.

The three time systems are described as follows:

TDB
    ("Barycentric Dynamical Time"); typically for cartesian, 
    osculating element, and close-approach tables. The uniform time 
    scale and independent variable of the planetary ephemeris 
    dynamical equations of motion.

TT
    ("Terrestrial (Dynamic) Time"), called TDT prior to 1991, used for 
    observer quantity tables. This is proper time as measured by 
    an Earth-bound observer and is directly related to atomic time, TAI. 
    TT periodically differs from TDB by, at most, 0.002 seconds.

UT
    is Universal Time This can mean one of two non-uniform time-scales 
    based on the rotation of the Earth. For this program, prior to 1962, 
    UT means UT1. After 1962, UT means UTC or "Coordinated Universal Time". 
    Future UTC leap-seconds are not known yet, so the closest known 
    leap-second correction is used over future time-spans. 
"""
