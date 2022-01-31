import cpuinfo
import pytest
"""
Functions to skip certain tests if not compatible with the computer 
architecture.
"""


def skip_m1(msg=None):
    """ Skip tests not compatible with the M1 Mac chip."""
    chip = cpuinfo.get_cpu_info().get('brand_raw')
    if 'm1' in chip.lower():
        reason = "does not work with M1 chips?"
        if msg is not None:
            reason = "{0} {1}".format(msg, reason)

        pytest.skip(reason)
