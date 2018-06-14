def chi2_telescope(your_event, model_1, parameters):
    chichi = 0.
    for telescope in your_event.telescopes:
        model = model_1.compute_the_microlensing_model(telescope, parameters)
        flux = telescope.lightcurve_flux[:,1]
        errflux = telescope.lightcurve_flux[:,2]
        residus = (flux - model[0])/errflux
        chichi += (residus ** 2).sum()
    return chichi

