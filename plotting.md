# Updating Plotting in MulensModel.

## Overview of the problem

There are three functions that plot data:

mulensdata.plot()
model.plot_data()
model.plot_residuals()

The desired behavior is that all three of these functions works the
same way except that they plot slightly different things. Therefore,
there should be a common set of arguments between them. As with the
event.plot_data() method, that common set should be passed with
**kwargs so they can be maintained in some centralized location. The
methods also need to take **kwargs for plotting.

Our special kwargs:
subtract_2450000, subtract_2460000
show_errorbars, show_bad

They all need to be able to access MulensData.plot_properties().

I think this means that we need a common "plot(t, x, **kwargs)"
function for all of them.

## Proposed solution

A new plot.py module that includes a plot.plot() method. 

*NOTE: After drafting these functions, there are clearly some
 inconsistencies/inefficiencies (e.g. where to handle
 separate_custom_kwargs). These should be reconsidered before
 implementing.*

### Outline of the plot.plot() method:

def plot(x, y, yerr=None, **kwargs):

   (custom_kwargs, plot_kwargs) = separate_custom_kwargs(**kwargs)
   x_sub = subtract(x, **custom_kwargs)
   plot_data(x_sub, y, yerr=yerr, **plot_kwargs)

def plot_data(x, y, yerr=None, **plot_kwargs):

    if yerr is not None:
       pl.errorbar(x, y, yerr=yerr, **plot_kwargs)
    else:
       pl.scatter(x, y, **plot_kwargs)

def set_plot_axes(**custom_kwargs):
    """
    xlabel
    ylabel
    handle subtract_245,246
    """
    pl.gca().invert_yaxis() - or maybe not as default? in case plotting in fluxes?
    
### Calling for model.plot_data(), model.plot_residuals(), mulensdata.plot():


#### model.py

Should plot_data have an option to plot in fluxes/

def plot_data(data_ref=None, show_errorbars=True, **kwargs):
    
    (f_source, f_blend) = get_ref_fluxes()
    for dataset in datasets:
    	flux = scale_flux(dataset, fluxes=(f_source, f_blend))
	mag = utils.flux2mag(flux)
	if show_errorbars:
  	    err_flux = scale_err_flux(dataset, fluxes=(f_source, f_blend))
	    err_mag = utils.flux_err2mag_err(flux, err_flux)
        else: 
	    err_flux = None

	new_kwargs = dataset.add_dataset_kwargs(**kwargs)
	plot.plot(dataset.times, mag, yerr=err_mag, **new_kwargs)

    plot.set_plot_axes(**kwargs)

def plot_residuals(data_ref=None, show_errorbars=True, **kwargs):
    """ Same as plot_data except with the appropriate residuals calculation"""
    
#### data.py

def plot(show_errorbars=True, phot_fmt=None, **kwargs):
    if phot_fmt is None:
       phot_fmt = dataset.phot_fmt

    yerr = None
    if phot_fmt == 'mag':
        y = dataset.mag	
	if show_errorbars:
	    yerr = dataset.err_mag
 
    elif phot_fmt == 'flux':
        y = dataset.flux
        if show_errorbars:
	    yerr = dataset.err_mag

    new_kwargs = dataset.add_dataset_kwargs(**kwargs)
    plot.plot(dataset.times, y, yerr=yerr, **new_kwargs)