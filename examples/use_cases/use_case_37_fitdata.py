"""
Use cases for the new FitData class
"""
import MulensModel as mm
import os

#############
# Define a basic PSPL fit
# Initial Model
t_0 = 2455379.571
u_0 = 0.5
t_E = 17.94

pspl_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Import data
file_name = os.path.join(mm.DATA_PATH, 'photometry_files',
                         'OB08092', 'phot_ob08092_O4.dat')
data = mm.MulensData(file_name=file_name)

# Freely fit the fluxes
free_fit = mm.FitData(model=pspl_model, dataset=data)
free_fit.fit_fluxes()

# Access different properties of that fit
print(free_fit.source_fluxes, free_fit.blend_flux)
#commented out since not implemented
'''
print(free_fit.chi2_per_datapoint)
print(free_fit.chi2)
'''

# since the rest hasn't been implemented, quitting early
quit()

# Constrain various aspects of the fit
fit_1 = mm.FitData(model=pspl_model, dataset=data, fix_source_flux=False, fix_blend_flux=0.)

fit_2 = mm.FitData(model=pspl_model, dataset=data, fix_blend_flux=0.1584)
# Could re-write to act on an OGLE DIA event with f_base forced to I=20.0

fit_3 = mm.FitData(model=pspl_model, dataset=data, fix_source_flux=0.1)
# Could re-write such that 0.1 --> known f_source value

#################
# Maybe add a version of the OB08092 data for which a binary source fit is plausible (i.e. including 2008 data)
# Define a basic binary source fit
binary_source_model = mm.Model(
    {'t_0_1': 2454727.39, 'u_0_1': 1.5502, 't_0_2': 2454541.208, 'u_0_2': 0.062, 't_E': 38.56})

# Freely fit the fluxes
free_1L2S_fit = mm.FitData(model=binary_source_model, dataset=data)

# Fix the flux ratio
free_1L2S_fit = mm.FitData(model=binary_source_model, dataset=data, fix_q_flux=0.629/38.56)