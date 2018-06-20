def chi2_telescope(your_event, model, parameters_list):
    parameters = model.compute_pyLIMA_parameters(parameters_list)
    chichi = 0.
    for telescope in your_event.telescopes:
        model_ = model.compute_the_microlensing_model(telescope, parameters)
        flux = telescope.lightcurve_flux[:,1]
        errflux = telescope.lightcurve_flux[:,2]
        residus = (flux - model_[0])/errflux
        chichi += (residus ** 2).sum()
    return chichi
