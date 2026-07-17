import os
import ctypes


try:
    import MulensModel.AdaptiveContouring as mm_ac
except Exception:
    _adaptive_contouring_wrapped = False
else:
    _adaptive_contouring_wrapped = True


def _try_load(path, name):
    """
    Try loading compiled C library. Input is *str* or *list* of *str*.
    """
    if isinstance(path, str):
        path = [path]

    for path_ in path:
        try:
            out = ctypes.cdll.LoadLibrary(path_)
        except OSError as error:
            print("WARNING - File not loaded:", path_)
            print("OSError:\n", error)
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


def _import_compiled_AdaptiveContouring():
    """try importing manually compiled AdaptiveContouring package"""
    ac = "AdaptiveContouring"
    adaptive_contour = _try_load(_get_path_2(ac, ac + "_wrapper.so"), ac)
    _adaptive_contouring_wrapped = (adaptive_contour is not None)
    if not _adaptive_contouring_wrapped:
        return (_adaptive_contouring_wrapped, None)
    adaptive_contour.Adaptive_Contouring_Linear.argtypes = (8 * [ctypes.c_double])
    adaptive_contour.Adaptive_Contouring_Linear.restype = ctypes.c_double
    return (_adaptive_contouring_wrapped, adaptive_contour.Adaptive_Contouring_Linear)


if _adaptive_contouring_wrapped:
    _adaptive_contouring_linear = mm_ac.Adaptive_Contouring_Linear
else:
    out = _import_compiled_AdaptiveContouring()
    _adaptive_contouring_wrapped = out[0]
    _adaptive_contouring_linear = out[1]
