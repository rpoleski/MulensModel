import numpy as np
from os.path import isfile

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy import __version__ as astropy_version

from MulensModel.utils import Utils

"""
All of the documentation Radek could find on Time reference frames and
scales (e.g. HJD_UTC vs. BJD_TDB) for astropy.Time and JPL Horizons
are copied to ../../documents/Horizons_notes.md.
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
        self._file_properties = {"file_name": file_name}
        self._read_input_file()

    def _read_input_file(self):
        """check the file type and then read it properly"""
        file_type = 'np.array'
        with open(self._file_properties['file_name'], 'r') as in_file:
            for line in in_file.readlines():
                if line[0:5] == '$$SOE':
                    file_type = 'Horizons'
                    break

        if not isfile(self._file_properties['file_name']):
            msg = 'Horizons files {:} does not exists.'
            message = msg.format(self._file_properties['file_name'])
            raise FileExistsError(message)

        if file_type == 'Horizons':
            self._read_horizons_file()
        else:
            (time, x, y, z) = np.loadtxt(
                self._file_properties['file_name'],
                usecols=(0, 1, 2, 3), unpack=True)
            self._time = time
            key = 'representation'
            if int(astropy_version[0]) >= 4:
                key = "representation_type"
            self._xyz = SkyCoord(x=x, y=y, z=z, **{key: 'cartesian'})

    def _get_start_end(self):
        """
        Find the start (self.start_ind) and end (self.stop_ind) points
        for the data in the Horizons file (i.e. the end of the header
        and beginning of the footer). Also counts the lines in the
        file. (self.line_count)
        """
        with open(self._file_properties['file_name'], 'r') as in_file:
            self._file_properties['line_count'] = 0
            for (i, line) in enumerate(in_file):
                # Check for "start data" string
                if line[0:5] == '$$SOE':
                    self._file_properties['start_ind'] = i

                # Check for "end data" string
                if line[0:5] == '$$EOE':
                    self._file_properties['stop_ind'] = i

                # Count total number of lines
                self._file_properties['line_count'] += 1

    def _read_horizons_file(self):
        """
        reads standard output from JPL Horizons
        """
        # Read in the file
        self._get_start_end()
        data = np.genfromtxt(
            self._file_properties['file_name'],
            dtype=[('date', 'S17'), ('ra_dec', 'S23'), ('distance', 'f8'),
                   ('foo', 'S23')],
            delimiter=[18, 29, 17, 24], autostrip=True,
            skip_header=self._file_properties['start_ind'] + 1,
            skip_footer=(self._file_properties['line_count'] -
                         self._file_properties['stop_ind']))

        # Fix time format
        for (i, date) in enumerate(data['date']):
            data['date'][i] = Utils.date_change(date)

        # Currently we assume HORIZONS works in UTC.
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
