import numpy as np

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from MulensModel.utils import Utils


month_3letter_to_2digit = {
        'Jan' : '01',
        'Feb' : '02',
        'Mar' : '03',
        'Apr' : '04',
        'May' : '05',
        'Jun' : '06',
        'Jul' : '07',
        'Aug' : '08',
        'Sep' : '09',
        'Oct' : '10',
        'Nov' : '11',
        'Dec' : '12'
        }

def date_change(text):
    """changes format: ' 2015-Oct-30 12:00' -> '2015-10-30 12:00' """
    n1 = text.find("-")
    n2 = text.rfind("-")
    if n1 == -1 or n1 == n2:
        raise ValueError("Can't run date_change() for {:}".format(text))
    return text[1:n1+1] + month_3letter_to_2digit[text[n1+1:n2]] + text[n2:]


class Horizons(object):
    def __init__(self, file_name=None):
        self._time = None
        self._xyz = None
        self._names = {}
        self._names['date'] = 'Date__(UT)__HR:MN'
        self._names['ra_dec'] = 'R.A._(ICRF/J2000.0)_DEC'
        self._names['distance'] = 'delta'
        apply_float = ['distance']
        apply_date_change = ['date']
        self.data_lists = {}
        for key in self._names:
            self.data_lists[key] = []
            
        self._read_horizons_file(file_name)
        
        for key in apply_float:
            for j in range(len(self.data_lists[key])):
                self.data_lists[key][j] = float(self.data_lists[key][j])
        for key in apply_date_change:
            for j in range(len(self.data_lists[key])):
                self.data_lists[key][j] = date_change(self.data_lists[key][j])

    def _read_horizons_file(self, file_name):
        """reads standard output from JPL Horizons"""
        with open(file_name) as in_file:
            mode = None # Other possible values are: 'header', 'header_done', 'main', and 'finished'.
            # This gives information in which part of the Horizons file we are currently in.
            for l in in_file.readlines():
                line = l[:-1]
                if len(line) == 0:
                    continue
                mode_save = mode
                if line[0] == '$':
                    if line == '$$SOE':
                        if mode == 'header_done':
                            mode = 'main'
                        else:
                            raise ValueError('error: {:}'.format(line))
                    elif line == '$$EOE':
                        if mode == 'main':
                            mode = 'finished'
                            break
                        else:
                            raise ValueError('error: {:}'.format(line))
                    else:
                        raise ValueError('unexpected input: {:}'.format(line))
                elif line[0] == '*':
                    if len(line) != 79:
                        if mode is None:
                            mode = 'header'
                        elif mode == 'header':
                            mode = 'header_done'
                    else:
                        if mode == 'main':
                            raise ValueError('error: {:}'.format(line))
                        
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
        """return HJD vector"""
        if self._time is None:
            self._time = Time(self.data_lists['date'], format='iso').jd - 2450000.
        return self._time

    @ property 
    def xyz(self):
        """return X,Y,Z positions"""
        if self._xyz is None:
            self._xyz = SkyCoord(self.data_lists['ra_dec'], distance=self.data_lists['distance'], unit=(u.hourangle, u.deg, u.au)).cartesian
        return self._xyz

