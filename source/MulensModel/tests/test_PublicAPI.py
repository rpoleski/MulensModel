import MulensModel as mm


def test_all_entries_are_resolvable():
    """
    Every name listed in MulensModel.__all__ must resolve to a real attribute
    of the package. Catches stale entries left behind after class renames
    (regression: 'BinaryLensPointSourceWithShearVBBLMagnification' was kept
    in __all__ after the VBBL -> VBM rename in commit 666bef44 even though
    no class with that name was exported, so `from MulensModel import *`
    crashed with AttributeError).
    """
    missing = [name for name in mm.__all__ if not hasattr(mm, name)]
    assert missing == [], (
        "MulensModel.__all__ references symbols that do not exist on the "
        "package: {}".format(missing))


def test_star_import_succeeds():
    """
    `from MulensModel import *` must not raise. CPython evaluates __all__
    during a star import and fails hard if any listed name is missing.
    """
    namespace = {}
    exec("from MulensModel import *", namespace)
    for name in mm.__all__:
        assert name in namespace, "{} missing from star-import namespace".format(name)
