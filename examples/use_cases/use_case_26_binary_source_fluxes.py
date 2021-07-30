"""
If a binary source is observed from more than one site, the two
components should have the same flux fraction if the band is the same.
"""
import MulensModel as mm
import scipy.optimize as op


raise NotImplementedError(
    'This use case has not been implemented. Needs fake data files.')


def chi2_fun(theta, event, parameters_to_fit):
    """
    for given event set attributes from parameters_to_fit (list of
    str) to values from theta list
    """
    for (theta_, parameter) in zip(theta, parameters_to_fit):
        if parameter[0] == 'q':
            setattr(event.fix_source_flux_ratio, parameter[2:], theta_)
        else:
            setattr(event.model.parameters, parameter, theta_)

    return event.get_chi2()


# Import Data
data_site1_band1 = mm.MulensData(file_name='DATA_FILE_1.dat', bandpass='I')
data_site1_band2 = mm.MulensData(file_name='DATA_FILE_2.dat', bandpass='V')
data_site2_band1 = mm.MulensData(file_name='DATA_FILE_3.dat', bandpass='I')
data_site2_band2 = mm.MulensData(file_name='DATA_FILE_4.dat', bandpass='V')
datasets = [
    data_site1_band1, data_site1_band2, data_site2_band1, data_site2_band2]

# Create a Model (totally made up)
model = mm.Model(
    parameters=dict(
        t_0_1=2458000., u_0_1=0.1, t_0_2=2458004., u_0_2=0.3, t_E=25.))

# Fit the model
event = mm.Event(datasets=datasets, model=model)

parameters_to_fit = ['t_0_1', 'u_0_1', 't_0_2', 'u_0_2', 't_E', 'q_I', 'q_V']
initial_guess = [
    model.parameters.t_0_1,
    model.parameters.u_0_1,
    model.parameters.t_0_2,
    model.parameters.u_0_2,
    model.parameters.t_E,
    0.3, 0.3]
# where q_I = f_source_2 / f_source_1 in I band (q_V is the equivalent
# in V band). Then, you would only do the linear fit for f_source_1
# and f_blend in each band.
result = op.minimize(
    chi2_fun, x0=initial_guess, args=(event, parameters_to_fit),
    method='Nelder-Mead')
print(result.x)
