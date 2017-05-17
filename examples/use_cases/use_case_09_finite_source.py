import MulensModel

#Delete? Redundant with UC18/19?

m = MulensModel.Model()
m.finite_source = True

m2 = MulensModel.Model()
m2.finite_source(method = 'Stokes', time_range = (7605., 7606.), check_range = True, tolerance = 0.0001) # time_range and check_range are optional


