#! /usr/bin/env python

import sys
import numpy as np
from astropy.time import Time

from MulensModel.model import Model
from MulensModel.mulensdata import MulensData


def test_model_PSPL_1():
    """tests basic evaluation of Paczynski model"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    data = MulensData()
    data._date_zeropoint = 2450000.
    times = np.array([t_0-2.5*t_E, t_0, t_0+t_E])
    data._time = Time(times+data._date_zeropoint, format="jd")
    model = Model(t_0=t_0)
    model.t_0 = t_0
    model.u_0 = u_0
    model.t_E = t_E
    model._datasets = [data]
    np.testing.assert_almost_equal(model.magnification, [np.array([1.028720763, 2.10290259, 1.26317278])], err_msg="PSPL model returns wrong values")

