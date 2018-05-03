from os import path

from .version import __version__
from MulensModel.event import *
from MulensModel.fit import *
from MulensModel.model import *
from MulensModel.modelparameters import *
from MulensModel.mulensdata import *
from MulensModel.binarylens import BinaryLens
from MulensModel.pointlens import PointLens

from MulensModel.mulensobjects import *


__all__ = ['mulensobjects', 'MODULE_PATH']

MODULE_PATH = path.abspath(__file__)
for i in range(3):
    MODULE_PATH = path.dirname(MODULE_PATH)
# We do the same in binarylens.py.
