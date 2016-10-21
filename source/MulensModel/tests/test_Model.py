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
    times = np.array([t_0-2.5*t_E, t_0, t_0+t_E])
    data = MulensData(data_list=[times, times*0., times*0.], date_fmt='jdprime')
    model = Model(t_0=t_0)
    model.u_0 = u_0
    model.t_E = t_E
    model.set_datasets([data])
    np.testing.assert_almost_equal(model.magnification, [np.array([1.028720763, 2.10290259, 1.26317278])], err_msg="PSPL model returns wrong values")

def test_model_init_1():
    """tests if basic parameters of Model.__init__() are properly passed"""
    t_0 = 5432.10987
    u_0 = 0.001
    t_E = 123.456
    rho = 0.0123
    m = Model(t_0=t_0, u_0=u_0, t_E=t_E, rho=rho)
    np.testing.assert_almost_equal(m.t_0, t_0, err_msg='t_0 not set properly')
    np.testing.assert_almost_equal(m.u_0, u_0, err_msg='u_0 not set properly')
    np.testing.assert_almost_equal(m.t_E, t_E, err_msg='t_E not set properly')
    np.testing.assert_almost_equal(m.rho, rho, err_msg='rho not set properly')

if __name__ == "__main__":
    test_model_init_1()

