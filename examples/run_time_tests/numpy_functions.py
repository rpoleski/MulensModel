from numpy import sqrt, square, reciprocal, sum, polyfit, ones, copy, array, linalg

def numpy_chi2_v1(time, obs_flux, obs_flux_err, t_0, u_0, t_E):
    u2 = ((time-t_0)/t_E)**2 + u_0**2
    magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
    fluxes = polyfit(magnification, obs_flux, 1, w=reciprocal(obs_flux_err))
    flux = fluxes[0] * magnification + fluxes[1]
    chi2 = sum(square((obs_flux - flux) / obs_flux_err))
    return chi2

def numpy_chi2_v2(time, obs_flux, obs_flux_err, t_0, u_0, t_E):
    u2 = ((time-t_0)/t_E)**2 + u_0**2
    magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
    x = ones(shape=(2, len(u2)))
    x[0] = magnification
    sigma_inverse = 1. / obs_flux_err
    y = obs_flux * sigma_inverse
    xT = (x * sigma_inverse).T
    results = linalg.lstsq(xT, y, rcond=-1)[0]
    flux = results[0] * magnification + results[1]
    chi2 = sum(square((obs_flux - flux) * sigma_inverse))
    return chi2

def get_magnification(time, t_0, u_0, t_E):
    u2 = ((time-t_0)/t_E)**2 + u_0**2
    magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
    return magnification

def get_fluxes(magnification, obs_flux, obs_flux_err):
    x = ones(shape=(2, len(magnification)))
    x[0] = magnification
    sigma_inverse = 1. / obs_flux_err
    y = obs_flux * sigma_inverse
    xT = (x * sigma_inverse).T
    results = linalg.lstsq(xT, y, rcond=-1)[0]
    return results

def numpy_chi2_v3(time, obs_flux, obs_flux_err, t_0, u_0, t_E):
    magnification = get_magnification(time, t_0, u_0, t_E)
    f = get_fluxes(magnification, obs_flux, obs_flux_err)
    flux = f[0] * magnification + f[1]
    chi2 = sum(square((obs_flux - flux) / obs_flux_err))
    return chi2
