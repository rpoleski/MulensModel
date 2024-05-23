from os import path

from MulensModel.binarylens import BinaryLens
from MulensModel.binarylenswithshear import BinaryLensWithShear
from MulensModel.caustics import Caustics
from MulensModel.causticspointwithshear import CausticsPointWithShear
from MulensModel.causticswithshear import CausticsWithShear
from MulensModel.coordinates import Coordinates
from MulensModel.event import Event
from MulensModel.fitdata import FitData
from MulensModel.horizons import Horizons
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.model import Model
from MulensModel.modelparameters import ModelParameters, which_parameters
from MulensModel.mulensdata import MulensData
from MulensModel.mulensobjects import Lens, Source, MulensSystem
from MulensModel import orbits
from MulensModel.pointlens import PointLens, get_pspl_magnification
from MulensModel.pointlenswithshear import PointLensWithShear
from MulensModel.pointlensfinitesource import PointLensFiniteSource
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.trajectory import Trajectory
from MulensModel.uniformcausticsampling import UniformCausticSampling
from MulensModel.utils import MAG_ZEROPOINT, Utils

from .version import __version__

__all__ = ['mulensobjects', 'MODULE_PATH', 'DATA_PATH', 'BinaryLens',
           'BinaryLensWithShear', 'Caustics', 'CausticsPointWithShear',
           'CausticsWithShear', 'Coordinates', 'Event', 'FitData', 'Horizons',
           'LimbDarkeningCoeffs', 'MagnificationCurve', 'Model',
           'ModelParameters', 'which_parameters', 'MulensData', 'Lens',
           'Source', 'MulensSystem', 'orbits', 'PointLens',
           'get_pspl_magnification', 'PointLensWithShear',
           'PointLensFiniteSource', 'SatelliteSkyCoord', 'Trajectory',
           'UniformCausticSampling', 'MAG_ZEROPOINT',  'Utils', '__version__']

MODULE_PATH = path.abspath(__file__)
for i in range(3):
    MODULE_PATH = path.dirname(MODULE_PATH)

path_1 = path.join(MODULE_PATH, 'data')
if path.isdir(path_1):
    DATA_PATH = path_1
else:
    DATA_PATH = path.join(path.dirname(__file__), 'data')
