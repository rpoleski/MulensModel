"""
Fit a binary source event. Allow the flux ratio to be freely fit for KMTC data,
but then constrained for other datasets in the same band.

For example, suppose there is a short-term anomaly that is only covered by KMTC
data. Then, for a binary source fit, the KMTC data constrain q_flux but the
other datasets do not. However, for a self-consistent fit, fits to KMTA and
KMTS data must still track the contribution from the second source because even
if it doesn't *appear* to contribute to the other datasets, it might. Possibly
this is not something you would ever want to do In Real Life, but the point is,
you could if you wanted to.

This use case is not functional. To make it functional, someone needs to track
down an event with appropriate data.
"""
import MulensModel as mm

raise NotImplementedError('Needs fake data.')

# define some fake data
files = ['phot.dat', 'KMTC_I.pysis', 'KMTA_I.pysis', 'KMTS_I.pysis',
         'KMTC_V.pysis', 'KMTA_V.pysis', 'KMTS_V.pysis']
bandpasses = ['I'] * 4 + ['V'] * 3
kwargs = {'phot_fmt': 'mag', 'usecols': range(3)}
datasets = [mm.MulensData(file_name=file_, bandpass=bandpass, **kwargs)
            for (file_, bandpass) in zip(files, bandpasses)]

# define the model
binary_source_model = mm.Model(
    {'t_0_1': 2459000.0, 'u_0_1': 1.5, 't_0_2': 2459007.0, 'u_0_2': 0.01,
     't_E': 30.})


# Create my own event class
class MyEvent(mm.Event):

    def fit_fluxes(self):
        """
        Allow the two source fluxes to be freely fit for some reference
        dataset, but then constrain the fluxes for all other datasets in
        the same bandpass.
        """
        self.fits = []
        kmtc_fits = {
            1: mm.FitData(model=self.model, dataset=self.datasets[1]),
            4: mm.FitData(model=self.model, dataset=self.datasets[4])}
        band = {'I': 1, 'V': 4}  # This simplies the code below.
        for (i, dataset) in enumerate(self.datasets):
            if i in kmtc_fits:
                fit = kmtc_fits[i]
            else:
                q_flux = kmtc_fits[band[dataset.bandpass]].source_flux_ratio
                fit = mm.FitData(model=self.model, dataset=dataset,
                                 fix_source_flux_ratio=q_flux)
            self.fits.append(fit)


# Fit the fluxes
event = MyEvent(model=binary_source_model, datasets=datasets)
print(event.chi2)
print(event.source_fluxes)
print(event.blend_fluxes)
