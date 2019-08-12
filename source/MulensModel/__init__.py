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


__all__ = ['mulensobjects', 'MODULE_PATH', 'DATA_PATH']

MODULE_PATH = path.abspath(__file__)
for i in range(3):
    MODULE_PATH = path.dirname(MODULE_PATH)

path_1 = path.join(MODULE_PATH, 'data')
if path.isdir(path_1):
    DATA_PATH = path_1
else:
    DATA_PATH = path.join(path.dirname(__file__), 'data')
