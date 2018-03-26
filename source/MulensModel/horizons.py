import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from MulensModel.utils import Utils

"""
All of the documentation Radek could find on Time reference frames and
scales (e.g. HJD_UTC vs. BJD_TDB) for astropy.Time and JPL Horizons
are copied at the end of this file.
"""


class Horizons(object):
    """
    An Object to read and hold the standard JPL Horizons output,
    i.e. satellite ephemerides.

    Arguments :
        file_name: *str*
            output from JPL Horizons file name

    For info on preparing JPL Horizons file for import, see instructions_.

    .. _instructions:
        https://github.com/rpoleski/MulensModel/blob/master/documents/Horizons_manual.md

    """

    def __init__(self, file_name):
        # initialize components
        self._time = None
        self._xyz = None

        # Read in the Horizons file
        self.file_properties = {"file_name": file_name}
        self._read_input_file()

    def _read_input_file(self):
        """check the file type and then read it properly"""
        file_type = 'np.array'
        with open(self.file_properties['file_name'], 'r') as in_file:
            for line in in_file.readlines():
                if line[0:5] == '$$SOE':
                    file_type = 'Horizons'
                    break

        if file_type == 'Horizons':
            self._read_horizons_file()
        else:
            (time, x, y, z) = np.loadtxt(self.file_properties['file_name'],
                usecols=(0, 1, 2, 3), unpack=True)
            self._time = time
            self._xyz = SkyCoord(x=x, y=y, z=z, representation='cartesian')

    def _get_start_end(self):
        """
        Find the start (self.start_ind) and end (self.stop_ind) points
        for the data in the Horizons file (i.e. the end of the header
        and beginning of the footer). Also counts the lines in the
        file. (self.line_count)
        """
        with open(self.file_properties['file_name'], 'r') as in_file:
            self.file_properties['line_count'] = 0
            for (i, line) in enumerate(in_file):
                # Check for "start data" string
                if line[0:5] == '$$SOE':
                    self.file_properties['start_ind'] = i

                # Check for "end data" string
                if line[0:5] == '$$EOE':
                    self.file_properties['stop_ind'] = i

                # Count total number of lines
                self.file_properties['line_count'] += 1

    def _read_horizons_file(self):
        """
        reads standard output from JPL Horizons
        """
        # Read in the file
        self._get_start_end()
        data = np.genfromtxt(
            self.file_properties['file_name'],
            dtype=[('date', 'S17'), ('ra_dec', 'S23'), ('distance', 'f8'),
                   ('foo', 'S23')],
            delimiter=[18, 29, 18, 24], autostrip=True,
            skip_header=self.file_properties['start_ind'] + 1,
            skip_footer=(self.file_properties['line_count'] -
                         self.file_properties['stop_ind']))

        # Fix time format
        for (i, date) in enumerate(data['date']):
            data['date'][i] = Utils.date_change(date)

        #  Currently we assume HORIZONS works in UTC.
        dates = [text.decode('UTF-8') for text in data['date']]
        self._time = Time(dates, format='iso', scale='utc').tdb.jd

        ra_dec = [text.decode('UTF-8') for text in data['ra_dec']]
        xyz = SkyCoord(
            ra_dec, distance=data['distance'], unit=(u.hourangle, u.deg, u.au))
        self._xyz = xyz.cartesian

    @property
    def time(self):
        """
        *np.ndarray*

        return times in TDB (reference frame depends on Horizons input)
        """
        return self._time

    @property
    def xyz(self):
        """
        *Astropy.CartesianRepresentation*

        return X,Y,Z positions based on RA, DEC and distance
        """
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
"""
