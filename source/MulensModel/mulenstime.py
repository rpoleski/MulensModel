from astropy.time import Time

class MulensTime(object):
    def __init__(self,time=0.,date_fmt="jd"):
        self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt)
        self._time = Time(time+self._date_zeropoint, format="jd")

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
        raise ValueError('Invalid value for date_fmt. Allowed values: "jd", "hjd", "jdprime", "hjdprime", "mjd"')

