from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u


class MulensTime(object):
    def __init__(self, time=0., date_fmt="jd"):
        self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt)
        earth_center = EarthLocation.from_geocentric(0., 0., 0., u.m)
        self._time = Time(time+self._date_zeropoint, format="jd",
                            location=earth_center)

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
        return self._time.jd

    @property
    def time(self):
        return self.jd - self._date_zeropoint

    def _get_date_zeropoint(self, date_fmt="jd"):
        """ Return the zeropoint of the date so it can be converted to
        the standard 245#### format."""
        if date_fmt == "jd" or date_fmt == "hjd":
            return 0.
        if date_fmt == "jdprime" or date_fmt == "hjdprime":
            return 2450000.
        if date_fmt == "mjd":
            return 2400000.5
        lst = '"jd", "hjd", "jdprime", "hjdprime", "mjd"'    
        raise ValueError('Invalid value for date_fmt. Allowed values: ' + lst)

