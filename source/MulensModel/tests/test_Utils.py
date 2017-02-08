import MulensModel.utils


def test_complex_fsum_1():
    z = [(0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1-1e+99j), (0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1+1e+99j), (0.1+0.1j)]
    assert MulensModel.utils.Utils.complex_fsum(z) == (1 + 0.8j)
