import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

import MulensModel as mm


SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "photometry_files", "OB08092", "phot_ob08092_O4.dat")


def test_model_coords():
    """
    Test Model.coords and different changes of format.
    """
    coords = SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg))
    model_1 = mm.Model({'t_0': 2450000, 'u_0': 0.1, 't_E': 100},
                       coords='18:00:00 -30:00:00')
    assert isinstance(model_1.coords, SkyCoord)
    assert model_1.coords.ra == coords.ra
    assert model_1.coords.dec == coords.dec
    assert model_1.coords.dec.deg == -30.00

    model_3 = mm.Model({'t_0': 2450000, 'u_0': 0.1, 't_E': 100})
    model_3.coords = '17:00:00 -27:32:14'
    assert model_3.coords.to_string('hmsdms') == '17h00m00s -27d32m14s'


def test_event_coords():
    """
    Test Event.coords and different changes of format.
    """
    coord_str_event = '15h30m00s +45d00m00s'
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.})

    event = mm.Event(datasets=data, model=model, coords='15:30:00 45:00:00')
    assert event.coords.to_string('hmsdms') == coord_str_event
    assert event.model.coords.to_string('hmsdms') == coord_str_event

    event_2 = mm.Event(coords='15:30:00 45:00:00')
    event_2.datasets = [data]
    assert event_2.coords.to_string('hmsdms') == coord_str_event

    event_3 = mm.Event(datasets=data, model=model, coords='15:30:00 45:00:00')
    event_3.model.coords = '5:10:15 20:25:30'
    new_coord_str = '05h10m15s +20d25m30s'
    assert event_3.model.coords.to_string('hmsdms') == new_coord_str


def check_event_coords(event, ra, dec):
    """
    For given Event instance event, check if .ra, .model.ra,
    .datasets[0].ra etc. are equal to ra and dec
    """
    np.testing.assert_almost_equal(event.coords.ra.value, ra)
    np.testing.assert_almost_equal(event.model.coords.ra.value, ra)
    np.testing.assert_almost_equal(event.datasets[0].coords.ra.value, ra)
    np.testing.assert_almost_equal(event.coords.dec.value, dec)
    np.testing.assert_almost_equal(event.model.coords.dec.value, dec)
    np.testing.assert_almost_equal(event.datasets[0].coords.dec.value, dec)


def test_event_coords_ra_dec_1():
    """checks passing coords from from Event to MulensData"""
    coords_str_1 = '03:00:00 +44:15:00'
    ra_1 = 45.
    dec_1 = 44.25

    ra_2_str = "35d"
    ra_2 = 35.
    dec_2_str = "+32:00:00"
    dec_2 = 32.

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.})
    event = mm.Event(datasets=data, model=model, coords=coords_str_1)

    check_event_coords(event, ra_1, dec_1)

    event.coords = '{0} {1}'.format(ra_2_str, dec_2_str)

    np.testing.assert_almost_equal(data.coords.ra.value, ra_2)
    np.testing.assert_almost_equal(data.coords.dec.value, dec_2)


def test_event_coords_ra_dec_2():
    """checks setting coords of Event after initialization"""
    ra_1_str = '01:00:00'
    dec_1_str = '+44:15:00'
    ra_1 = 15.
    dec_1 = 44.25

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.})
    event = mm.Event(datasets=data, model=model)
    event.coords = '{0} {1}'.format(ra_1_str, dec_1_str)
    check_event_coords(event, ra_1, dec_1)


def test_v_Earth_projected():
    """checks function that calculates Earth projected velocity"""
    # Yee et al. 2015, ob140939:
    coords = mm.Coordinates("17:47:12.25 -21:22:58.7")
    np.testing.assert_almost_equal(
        [-0.5, 28.9],
        coords.v_Earth_projected(2456836.06), decimal=1)

    # Batista et al. 2011, mb09387:
    coords = mm.Coordinates("17:53:50.79 -33:59:25")
    np.testing.assert_almost_equal(
        [-3.60, 22.95],
        coords.v_Earth_projected(2455042.34,), decimal=2)
