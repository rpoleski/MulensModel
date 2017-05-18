
import MulensModel


# Prepare a list of 3 MulensData objects. Two of them have bandpass specified:
data_1 = MulensModel.MulensData(file_name="ogle_data_phot.dat", bandpass='I')
data_2 = MulensModel.MulensData(file_name="kmt_V_phot.dat", bandpass='V')

# and one has only relative weights of limb darkening coefficients:
data_3 = MulensModel.MulensData(file_name="moa_phot.dat")
moa_LD_weights = {'I': 2., 'V': 1.}
data_3.set_limb_darkening_weights(moa_LD_weights) 

# Initialize instance of Model.
model = MulensModel.Model(n_components=1)
model.set_datasets([data_1, data_2, data_3])

# Set coefficients - you can use either u (most common outside microlensing)
# or gamma (fixed total flux) coefficients. 
model.set_limb_coef_gamma('I', 0.4555)
model.set_limb_coef_u('V', 0.5555)

# Print LD coefficients:
print("bandpass   u   gamma:")
print("---------------------")
for band in model.bandpasses:    
    print("{:} {:.4f} {:.4f}".format(band, model.limb_coef_u(band), 
                                        model.limb_coef_gamma(band)))

#Add implementation for setting method to calculate FS effects for point lens
