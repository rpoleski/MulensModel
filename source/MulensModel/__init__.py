from os import path

from MulensModel.binarylens import BinaryLensPointSourceWM95Magnification,\
    BinaryLensPointSourceVBBLMagnification, BinaryLensPointSourceMagnification, \
    BinaryLensQuadrupoleMagnification, BinaryLensHexadecapoleMagnification, \
    BinaryLensVBBLMagnification, BinaryLensAdaptiveContouringMagnification
from MulensModel.binarylenswithshear import \
    BinaryLensPointSourceWithShearWM95Magnification, \
    BinaryLensPointSourceWithShearVBBLMagnification
from MulensModel.causticsbinary import CausticsBinary
from MulensModel.causticspointwithshear import CausticsPointWithShear
from MulensModel.causticsbinarywithshear import CausticsBinaryWithShear
from MulensModel.coordinates import Coordinates
from MulensModel.event import Event
from MulensModel.fitdata import FitData
from MulensModel.horizons import Horizons
from MulensModel.limbdarkeningcoeffs import LimbDarkeningCoeffs
from MulensModel.magnificationcurve import MagnificationCurve
from MulensModel.model import Model
from MulensModel.modelparameters import ModelParameters
from MulensModel.mulensdata import MulensData
from MulensModel.mulensobjects import Lens, Source, MulensSystem
from MulensModel import orbits
from MulensModel.pointlens import PointSourcePointLensMagnification, \
    FiniteSourceUniformGould94Magnification, FiniteSourceLDYoo04Magnification
from MulensModel.pointlenswithshear import PointSourcePointLensWithShearMagnification
from MulensModel.b0b1utils import B0B1Utils
from MulensModel.elliputils import EllipUtils
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.trajectory import Trajectory
from MulensModel.uniformcausticsampling import UniformCausticSampling
from MulensModel.utils import MAG_ZEROPOINT, Utils

from .version import __version__

__all__ = [
    'BinaryLensPointSourceWM95Magnification', 'BinaryLensPointSourceVBBLMagnification',
    'BinaryLensQuadrupoleMagnification', 'BinaryLensHexadecapoleMagnification', 'BinaryLensVBBLMagnification',
    'BinaryLensAdaptiveContouringMagnification', 'BinaryLensPointSourceWithShearWM95Magnification',
    'BinaryLensPointSourceWithShearVBBLMagnification', 'CausticsBinary', 'CausticsPointWithShear',
    'CausticsBinaryWithShear', 'Coordinates', 'Event', 'FitData', 'Horizons', 'LimbDarkeningCoeffs',
    'MagnificationCurve', 'Model', 'ModelParameters', 'MulensData', 'Lens', 'Source', 'MulensSystem', 'orbits',
    'PointSourcePointLensMagnification', 'FiniteSourceUniformGould94Magnification',
    'FiniteSourceLDYoo04Magnification', 'PointSourcePointLensWithShearMagnification', 'B0B1Utils', 'EllipUtils',
    'SatelliteSkyCoord', 'Trajectory', 'UniformCausticSampling', 'MAG_ZEROPOINT', 'Utils', '__version__']

# Set MODULE_PATH to the module directory
MODULE_PATH = path.dirname(path.abspath(__file__))

# Set DATA_PATH to the data directory within the package
DATA_PATH = path.join(path.dirname(__file__), 'data')
