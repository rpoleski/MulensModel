import MulensModel as mm


raise NotImplementedError('input data for this use case do not exist.')

# Prepare a list of 3 MulensData objects. Two of them have bandpass specified:
data_1 = mm.MulensData(file_name="ogle_data_phot.dat", bandpass='I')
data_2 = mm.MulensData(file_name="kmt_V_phot.dat", bandpass='V')

# and one has only relative weights of limb darkening coefficients:
data_3 = mm.MulensData(file_name="moa_phot.dat")
moa_LD_weights = {'I': 2., 'V': 1.}
data_3.set_limb_darkening_weights(moa_LD_weights)

# Initialize instance of Model.
model = mm.Model({'t_0': 2457600., 'u_0': 0.01, 't_E': 3.14})
model.set_datasets([data_1, data_2, data_3])

# Set finite source method (no limb darkening).
methods = [2455746., 'finite_source_uniform_Gould94', 2455746.6]
model.set_magnification_methods(methods)

# Set coefficients - you can use either u (most common outside microlensing)
# or gamma (fixed total flux) coefficients.
model.set_limb_coeff_gamma('I', 0.4555)
model.set_limb_coeff_u('V', 0.5555)

# Print LD coefficients:
print("bandpass   u   gamma:")
print("---------------------")
for band in model.bandpasses:
    print("{:} {:.4f} {:.4f}".format(band, model.get_limb_coeff_u(band),
                                     model.get_limb_coeff_gamma(band)))

# Add implementation for setting method to calculate FS effects for point lens
