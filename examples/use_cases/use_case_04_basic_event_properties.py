from astropy import units as u

import MulensModel


event = MulensModel.Event()
event.model = MulensModel.Model(t_0=7600., u_0=0.01, t_E=34.*u.day)

event.model.ra = 270.12345 * u.deg
event.model.dec = -27.6789 * u.deg

event.datasets.append(
    MulensModel.MulensData(np.array(hjd), np.array(m), np.array(e)))    
# where t, m, e are previously read/defined numpy arrays

event.datasets.append(
    MulensModel.MulensData(file_name='XXX.dat', format='OGLE', bandpass='I')) 
# possible formats: pysis, MOA, KMT, normal, OR pass a dictionary as below

event.datasets.append(
    MulensModel.MulensData(
        file_name='other_data.dat', 
        format=dict('hjd': 3, 'err_mag': 8, 'mag': 7))

