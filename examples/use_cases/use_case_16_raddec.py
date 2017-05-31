import MulensModel 
from astropy.coordinates import SkyCoord
from astropy import units as u
"""
Use cases for passing RA, DEC to MulensDAta, Model, and Event.

Based on OGLE-2015-BLG-0448 from Poleski et al. 2016 ApJ 823, 63
"""

#Specifying coordinates to calculate HJD from JD
data_1 = MulensModel.MulensData(
    file_name='../../data/ob150448_OGLE_ref_v1.dat', 
    coords='18:10:14.38 -31:45:09.4')

data_2 = MulensModel.MulensData(
    file_name='../../data/ob150448_OGLE_ref_v1.dat', ra='18:10:14.38', 
    dec='-31:45:09.4')

coords = SkyCoord('18:10:14.38 -31:45:09.4', unit=(u.hourangle, u.deg))
data_3 = MulensModel.MulensData(
    file_name='../../data/ob150448_OGLE_ref_v1.dat', coords=coords)

#Specifiying coordinates to calculate a model with parallax
t_0 = 2457213.146
u_0 = -0.0874
t_E = 61.02
pi_E_N = 0.1142
pi_E_E = -0.1088

ground_model = MulensModel.Model()
ground_model.set_parameters(t_0=t_0, u_0=u_0, t_E=t_E, pi_E=[pi_E_N, pi_E_E],
                 coords='18:10:14.38 -31:45:09.4')
space_model = MulensModel.Model(
    t_0=t_0, u_0=u_0, t_E=t_E, pi_E=[pi_E_N, pi_E_E], 
    ra='18:10:14.38', dec='-31:45:09.4', satellite='Spitzer', 
    ephemrides_file='Spitzer_ephemeris_01.dat')

#Access Galactic and ecliptic coordinates:
print(ground_model.galactic_l)
print(ground_model.galactic_b)
print(ground_model.ecliptic_lon)
print(ground_model.ecliptic_lat)

ground_model.plot_magnification(label='ground')
space_model.plot_magnification(label='space')
pl.legend()

#Sepcifying coordinates for an event
ground_data = MulensModel.MulensData(
    file_name='../../data/ob150448_OGLE_ref_v1.dat')
space_data = MulensModel.MulensData(
    file_name='../../data/ob150448_Spitzer_ref_v1.dat', 
    satellite='Spitzer', ephemrides_file='Spitzer_ephemeris_01.dat')

model_params = MulensModel.ModelParameters(
    t_0=t_0, u_0=u_0, t_E=t_E, pi_E_N=pi_E_N, pi_E_E=pi_E_E)
event = MulensModel.Event(datasets=[my_data], 
                          model=MulensModel.Model(parameters=model_params), 
                          coords='18:10:14.38 -31:45:09.4')
#Given these three different cases, it will be possible to specify
#conflicting sets of coordinates at different stages of model
#definition. Do we care?

pl.show()
