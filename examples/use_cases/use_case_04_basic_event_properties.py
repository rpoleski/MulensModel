import MulensModel
from astropy import units as u

e = MulensModel.Event()
e.model = MulensModel.Model(...)    # define in one of the ways from previous use cases

e.ra = 270.12345 * u.deg
e.dec = -27.6789 * u.deg

e.datasets.append(np.array(hjd), np.array(m), np.array(e))    # where t, m, e are previously read/defined numpy arrays

e.datasets.append(file_name='XXX.dat', format='OGLE', bandpass='I') # possible formats: pysis, MOA, KMT, normal, OR pass a dictionary as below

e.datasets.append(file_name='other_data.dat', format=dict('hjd': 3, 'err_mag': 8, 'mag': 7)

