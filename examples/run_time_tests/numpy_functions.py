from numpy import sqrt, square, reciprocal, sum, polyfit

def numpy_chi2_v1(time, obs_flux, obs_flux_err, t_0, u_0, t_E):
    u2 = ((time-t_0)/t_E)**2 + u_0**2
    magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
    fluxes = polyfit(magnification, obs_flux, 1, w=reciprocal(obs_flux_err))
    flux = fluxes[0] * magnification + fluxes[1]
    chi2 = sum(square((obs_flux - flux) / obs_flux_err))
    return chi2
