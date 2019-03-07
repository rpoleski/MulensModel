import MulensModel


s = "0.8 0.1 0.01 0.01 0.01 0.0 0.001"
f = [float(ss) for ss in s.split()]
print(MulensModel.binarylens._vbbl_binary_mag_dark(*f))

