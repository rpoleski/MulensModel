from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

from MulensModel.model import Model
from MulensModel.mulensdata import MulensData

def test_model_coords():
    coords = SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg))
    model_1 = Model(coords='18:00:00 -30:00:00')
    assert isinstance(model_1.coords, SkyCoord)
    assert model_1.coords.ra == coords.ra
    assert model_1.coords.dec == coords.dec
    assert  model_1.dec.deg == -30.00

    ra_2 = '17:00:00'
    dec_2 = '40:03:01'
    coords_2 = SkyCoord(
        '{0} {1}'.format(ra_2, dec_2), unit=(u.hourangle, u.deg))
    model_2 = Model()
    model_2.ra = ra_2
    model_2.dec = dec_2
    assert model_2.coords.ra == coords_2.ra
    assert model_2.coords.dec == coords_2.dec
    assert model_2.coords.to_string('hmsdms') == '17h00m00s +40d03m01s'

    model_3 = Model()
    model_3.coords = '17:00:00 -27:32:14'
    assert model_3.coords.to_string('hmsdms') == '17h00m00s -27d32m14s'
