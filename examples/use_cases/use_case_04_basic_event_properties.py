from astropy import units as u

import MulensModel


e = MulensModel.Event()
e.parameters(t_0=7603.1, u_0=0.23, t_E=45*u.day, rho=0.001, alpha=130.23*u.deg, s=1.3, q=0.3)

e.ra = 270.12345 * u.deg
e.dec = -27.6789 * u.deg

e.datasets.append(MulensModel.MulensData(np.array(hjd), np.array(m), np.array(e)))    # where t, m, e are previously read/defined numpy arrays

e.datasets.append(MulensModel.MulensData(file_name='XXX.dat', format='OGLE', bandpass='I')) # possible formats: pysis, MOA, KMT, normal, OR pass a dictionary as below

e.datasets.append(MulensModel.MulensData(file_name='other_data.dat', format=dict('hjd': 3, 'err_mag': 8, 'mag': 7))

