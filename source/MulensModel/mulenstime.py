from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils


zeropoints = {}
zeropoints["jd"] = 0.
zeropoints["rjd"] = 2400000.
zeropoints["mjd"] = 2400000.5
zeropoints["jdprime"] = 2450000.
zeropoints["hjd"] = zeropoints["jd"]
zeropoints["rhjd"] = zeropoints["rjd"]
zeropoints["hjdprime"] = zeropoints["jdprime"]


class MulensTime(object):
    def __init__(self, time=0., date_fmt="jd", coords=None):
        self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt.lower())
        if 'hjd' in date_fmt.lower():
            self._time_type = 'hjd'
        elif 'jd' in date_fmt.lower():
            self._time_type = 'jd'
        else:
            raise ValueError('Unrecognized time type (HJD/JD)')
        earth_center = EarthLocation.from_geocentric(0., 0., 0., u.m)
        self._time = Time(time+self._date_zeropoint, format="jd",
                            location=earth_center)
        if hasattr(time, '__iter__'):
            if any(self._time.jd < 2448600): # This corresponds to Dec 1991.
                raise ValueError('The time set is too early: {:}'.format(self._time.jd))
            elif any(self._time.jd > 2470000): # This corresponds to Jul 2050.
                raise ValueError('The time set is too late: {:}'.format(self._time.jd))
        else:
            if self._time.jd < 2448600 and time != 0.:
                raise ValueError('The time set is too early: {:}'.format(self._time.jd))
            elif self._time.jd > 2470000 and time != 0.:
                raise ValueError('The time set is too late: {:}'.format(self._time.jd))
        self._time_corr = None
        self._target = None
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._target = coords
            else:
                msg = 'unsupported format of coords parameter in MulensData()'
                raise ValueError(msg)        

    @property
    def astropy_time(self):
        """return astropy.Time object"""
        return self._time

    @property
    def zeropoint(self):
        """return the zeropoint of time"""
        return self._date_zeropoint

    @property
    def jd(self):
        """return full Julian Date"""
        if self._time_type == 'jd':
            return self._time.jd
        else:
            return (self._time - self._time_correction).jd
        
    @property
    def hjd(self):
        """return full Heliocentric Julian Date"""
        if self._time_type == 'hjd':
            return self._time.jd
        else:
            return (self._time + self._time_correction).jd        

    @property
    def _time_correction(self):
        '''time correction: HJD = JD + corr'''
        if self._time_corr is None:
            if self._target is None:
                msg1 = 'Event coordinates in MulensTime not set.\n'
                msg2 = "They're needed to calculate JD-HJD correction."
                raise ValueError(msg1 + msg2)
            if not Utils.astropy_version_check('1.2'):
                raise ValueError('Astropy 1.2 required for Time.light_travel_time()')
            self._time_corr = self._time.light_travel_time(self._target, 
                            'heliocentric')
        return self._time_corr
        
    @property
    def time(self):
        return self.jd - self._date_zeropoint

    def _get_date_zeropoint(self, date_fmt="jd"):
        """ Return the zeropoint of the date so it can be converted to
        the standard 245#### format."""
        if date_fmt not in zeropoints.keys():
            lst = '"' + '", "'.join(list(zeropoints.keys())) + '"'
            raise ValueError('Invalid value for date_fmt. Allowed values: ' + lst)
        return zeropoints[date_fmt]
