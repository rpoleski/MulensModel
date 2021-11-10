import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

import MulensModel as mm


def test_coords_format():
    """
    Simple test of __init__ and calculation of .x and .y.
    It checks if coordinates are passed properly in different formats.
    """
    times = np.linspace(2456789.0123, 2468013.579)

    parameters = {'t_0': 2462401., 'u_0': 0.1, 't_E': 2000.}
    params = mm.ModelParameters(parameters)

    coords_txt = "18:18:18.18 -30:30:30.30"
    coords = [None, coords_txt, mm.Coordinates(coords_txt),
              SkyCoord(coords_txt, unit=(u.hourangle, u.deg))]

    trajecotories = [mm.Trajectory(times, params, coords=c) for c in coords]
    for trajectory in trajecotories:
        assert np.all(trajectory.x == trajecotories[0].x)
        assert np.all(trajectory.y == trajecotories[0].y)

    parameters['pi_E_E'] = 0.1
    parameters['pi_E_N'] = -0.15
    params = mm.ModelParameters(parameters)

    coords = coords[1:]
    p = {'earth_orbital': True, 'satellite': False, 'topocentric': False}
    kwargs = {'times': times, 'parameters': params, 'parallax': p}
    trajecotories = [mm.Trajectory(coords=c, **kwargs) for c in coords]
    for trajectory in trajecotories:
        assert np.all(trajectory.x == trajecotories[0].x)
        assert np.all(trajectory.y == trajecotories[0].y)
