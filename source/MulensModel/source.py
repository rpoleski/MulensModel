from astropy import units as u


class Source(object):
    """
    Physical properties of a source (background) star.
    """
    def __init__(distance=None, angular_size=None,
                 limb_darkening=None):
        self.distance=distance
        self.angular_size=angular_size
        if limb_darkening is None:
            self.limb_darkening = LimbDarkeningCoeffs()
