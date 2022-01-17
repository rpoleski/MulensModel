import os
import ctypes
import numpy as np

try:
    import MulensModel.VBBL as mm_vbbl
except Exception:
    _vbbl_wrapped = False
else:
    _vbbl_wrapped = True
try:
    import MulensModel.AdaptiveContouring as mm_ac
except Exception:
    _adaptive_contouring_wrapped = False
else:
    _adaptive_contouring_wrapped = True


def _try_load(path, name):
    """
    Try loading compiled C library.
    Input is *str* or *list* of *str*.
    """
    if isinstance(path, str):
        path = [path]
    for path_ in path:
        try:
            out = ctypes.cdll.LoadLibrary(path_)
        except OSError:
            print("WARNING - File not loaded:", path_)
            print("Everything should work except:", name)
            pass
        else:
            return out
    return None


def _get_path_2(name_1, name_2):
    """convenience function"""
    module_path = os.path.abspath(__file__)
    for i in range(3):
        module_path = os.path.dirname(module_path)
    return os.path.join(module_path, 'source', name_1, name_2)


def _import_compiled_VBBL():
    """try importing manually compiled VBBL package"""
    vbbl = _try_load(
        _get_path_2('VBBL', "VBBinaryLensingLibrary_wrapper.so"), "VBBL")
    _vbbl_wrapped = (vbbl is not None)
    if not _vbbl_wrapped:
        return (_vbbl_wrapped, None, None, None, None)

    vbbl.VBBinaryLensing_BinaryMagDark.argtypes = 7 * [ctypes.c_double]
    vbbl.VBBinaryLensing_BinaryMagDark.restype = ctypes.c_double

    vbbl.VBBinaryLensing_BinaryMag0.argtypes = 7 * [ctypes.c_double]
    vbbl.VBBinaryLensing_BinaryMag0.restype = ctypes.c_double

    vbbl.VBBL_SG12_5.argtypes = 12 * [ctypes.c_double]
    vbbl.VBBL_SG12_5.restype = np.ctypeslib.ndpointer(
        dtype=ctypes.c_double, shape=(10,))

    vbbl.VBBL_SG12_9.argtypes = 20 * [ctypes.c_double]
    vbbl.VBBL_SG12_9.restype = np.ctypeslib.ndpointer(
        dtype=ctypes.c_double, shape=(18,))

    return (_vbbl_wrapped,
            vbbl.VBBinaryLensing_BinaryMagDark, vbbl.VBBL_SG12_5,
            vbbl.VBBinaryLensing_BinaryMag0, vbbl.VBBL_SG12_9)


def _import_compiled_AdaptiveContouring():
    """try importing manually compiled AdaptiveContouring package"""
    ac = "AdaptiveContouring"
    adaptive_contour = _try_load(_get_path_2(ac, ac + "_wrapper.so"), ac)
    _adaptive_contouring_wrapped = (adaptive_contour is not None)
    if not _adaptive_contouring_wrapped:
        return (_adaptive_contouring_wrapped, None)
    adaptive_contour.Adaptive_Contouring_Linear.argtypes = (
        8 * [ctypes.c_double])
    adaptive_contour.Adaptive_Contouring_Linear.restype = ctypes.c_double
    return (_adaptive_contouring_wrapped,
            adaptive_contour.Adaptive_Contouring_Linear)


# Check import and try manually compiled versions.
if _vbbl_wrapped:
    _vbbl_binary_mag_dark = mm_vbbl.VBBinaryLensing_BinaryMagDark
    _vbbl_binary_mag_0 = mm_vbbl.VBBinaryLensing_BinaryMag0
    _vbbl_SG12_5 = mm_vbbl.VBBL_SG12_5
    _vbbl_SG12_9 = mm_vbbl.VBBL_SG12_9
else:
    out = _import_compiled_VBBL()
    _vbbl_wrapped = out[0]
    _vbbl_binary_mag_dark = out[1]
    _vbbl_SG12_5 = out[2]
    _vbbl_binary_mag_0 = out[3]
    _vbbl_SG12_9 = out[4]

if not _vbbl_wrapped:
    _solver = 'numpy'
else:
    _solver = 'Skowron_and_Gould_12'
if _adaptive_contouring_wrapped:
    _adaptive_contouring_linear = mm_ac.Adaptive_Contouring_Linear
else:
    out = _import_compiled_AdaptiveContouring()
    _adaptive_contouring_wrapped = out[0]
    _adaptive_contouring_linear = out[1]
