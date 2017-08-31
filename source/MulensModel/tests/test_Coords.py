import sys, os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

import MulensModel
from MulensModel.model import Model
from MulensModel.mulensdata import MulensData
from MulensModel.event import Event


SAMPLE_FILE_01 = os.path.join(MulensModel.MODULE_PATH, 
                                    "data", "phot_ob08092_O4.dat")


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

def test_data_coords():
    coords = SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg))
    data_1 = MulensData(
        file_name=SAMPLE_FILE_01,
        coords='18:00:00 -30:00:00')
    assert isinstance(data_1.coords, SkyCoord)
    assert data_1.coords.ra == coords.ra
    assert data_1.coords.dec == coords.dec
    assert data_1.dec.deg == -30.00

    ra_2 = '17:00:00'
    dec_2 = '40:03:01'
    coords_2 = SkyCoord(
        '{0} {1}'.format(ra_2, dec_2), unit=(u.hourangle, u.deg))
    data_2 = MulensData(file_name=SAMPLE_FILE_01)
    data_2.ra = ra_2
    data_2.dec = dec_2
    assert data_2.coords.ra == coords_2.ra
    assert data_2.coords.dec == coords_2.dec
    assert data_2.coords.to_string('hmsdms') == '17h00m00s +40d03m01s'

    data_3 = MulensData(file_name=SAMPLE_FILE_01)
    data_3.coords = '17:00:00 -27:32:14'
    assert data_3.coords.to_string('hmsdms') == '17h00m00s -27d32m14s'

def test_event_coords():
    coord_str_event = '15h30m00s +45d00m00s'
    data = MulensData(
        file_name=SAMPLE_FILE_01,
        coords='00:00:15 -75:30:15')
    model = Model()

    event = Event(datasets=data, model=model, coords='15:30:00 45:00:00')
    assert event.coords.to_string('hmsdms') == coord_str_event
    assert event.model.coords.to_string('hmsdms') == coord_str_event
    assert event.datasets[0].coords.to_string('hmsdms') == coord_str_event

    coord_str_data = '00h00m15s -75d30m15s'
    data.coords = '00:00:15 -75:30:15'
    event_2 = Event(coords='15:30:00 45:00:00')
    event_2.datasets = [data]
    assert event_2.coords.to_string('hmsdms') == coord_str_data
    assert event_2.datasets[0].coords.to_string('hmsdms') == coord_str_data

    event_3 = Event(datasets=data, model=model, coords='15:30:00 45:00:00')
    event_3.model.coords = '5:10:15 20:25:30'
    new_coord_str = '05h10m15s +20d25m30s'
    assert event_3.model.coords.to_string('hmsdms') == new_coord_str
    #assert event_3.coords.to_string('hmsdms') == new_coord_str
    #assert event_3.datasets[0].coords.to_string('hmsdms') == new_coord_str
    #I don't think this worked previously, and I don't think it should
    #be allowed. - JCY

def check_event_coords(event, ra, dec):
    """For given Event instance event, check if .ra, .model.ra, 
    .datasets[0].ra etc. are equal to ra and dec"""
    np.testing.assert_almost_equal(event.ra.value, ra)
    np.testing.assert_almost_equal(event.model.ra.value, ra)
    np.testing.assert_almost_equal(event.datasets[0].ra.value, ra)
    np.testing.assert_almost_equal(event.dec.value, dec)
    np.testing.assert_almost_equal(event.model.dec.value, dec)
    np.testing.assert_almost_equal(event.datasets[0].dec.value, dec)


def test_event_coords_ra_dec_1():
    coords_str_1 = '03:00:00 +44:15:00'
    ra_1 = 45.
    dec_1 = 44.25

    ra_2_str = "35d"
    ra_2 = 35.
    dec_2_str = "+32:00:00"
    dec_2 = 32.

    data = MulensData(file_name=SAMPLE_FILE_01)
    model = Model()
    event = Event(datasets=data, model=model, coords=coords_str_1)

    # Assertations start here.
    check_event_coords(event, ra_1, dec_1)
    
    event.ra = ra_2_str
    check_event_coords(event, ra_2, dec_1)

    event.dec = dec_2_str
    check_event_coords(event, ra_2, dec_2)

    np.testing.assert_almost_equal(data.ra.value, ra_2)
    np.testing.assert_almost_equal(data.dec.value, dec_2)

def test_event_coords_ra_dec_2():
    ra_1_str = '01:00:00'
    dec_1_str = '+44:15:00'
    ra_1 = 15.
    dec_1 = 44.25

    data = MulensData(file_name=SAMPLE_FILE_01)
    model = Model()
    event = Event(datasets=data, model=model)
    event.ra = ra_1_str
    event.dec = dec_1_str

    check_event_coords(event, ra_1, dec_1)

