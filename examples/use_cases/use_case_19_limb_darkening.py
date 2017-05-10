
import MulensModel


# Prepare a list of 3 MulensData objects. Two of them have bandpass specified, 
# and one has only relative weights of limb darkening coefficients.
data_1 = MulensModel.MulensData(file_name="ogle_data_phot.dat", format='OGLE', bandpass='I')
data_2 = MulensModel.MulensData(file_name="moa_phot.dat", format="MOA")
moa_LD_weights = {'I': 2., 'V': 1.}
data_2.set_limb_darkening_weights(moa_LD_weights)
data_3 = MulensModel.MulensData(file_name="kmt_V_phot.dat", format="KMT", banpass='V')
datasets = [data_1, data_2, data_3]

# Initialize instance of Model.
model = MulensModel.Model(n_components=1)
model.set_datasets(datasets)

# Set coefficients - you can use either u (most common outside microlensing)
# or gamma (fixed total flux) coefficients. 
model.set_limb_coef_gamma('I', 0.4555)
model.set_limb_coef_u('V', 0.5555)

# Print bandpasses and LD coefs:
print("Dataset   bandpass   u   gamma")
for (i, data) in enumerate(datasets):
    print("{:} {:.4f} {:.4f}".format(i, data.bandpass, data.limb_coef_u, data.limb_coef_gamma))

