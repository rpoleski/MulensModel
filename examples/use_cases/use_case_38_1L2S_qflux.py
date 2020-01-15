"""
Fit a binary source event. Allow the flux ratio to be freely fit for KMTC data,
but then constrained for other datasets in the same band.

This use case is not functional. To make it functional, someone needs to track
down an event with appropriate data.
"""
import MulensModel as mm

# define some fake data
files = ['phot.dat', 'KMTC_I.pysis', 'KMTA_I.pysis', 'KMTS_I.pysis',
         'KMTC_V.pysis', 'KMTA_V.pysis', 'KMTS_V.pysis']
bandpasses = ['I'] * 4 + ['V'] * 3
kwargs = {'phot_fmt': 'mag', 'usecols': range(3)}
datasets = [mm.MulensData(file_name=file_, bandpass=bandpass, **kwargs)
            for (file_, bandpass) in zip(files, bandpasses)]

# define the model
binary_source_model = mm.Model(
    {'t_0_1': 2459000.0, 'u_0_1': 1.5, 't_0_2': 2459007.0, 'u_0_2': 0.01, 't_E': 30.})

# Create my own event class
class MyEvent(mm.Event):

    def fit_fluxes(self):
        """
        Allow the two source fluxes to be freely fit for some reference
        dataset, but then constrain the fluxes for all other datasets in
        the same bandpass.
        """
        self.fits = []
        kmtc_i_fit = mm.FitData(model=self.model, dataset=self.datasets[1])
        kmtc_v_fit = mm.FitData(model=self.model, dataset=self.datasets[4])
        for i, dataset in self.datasets:
            if i == 1:
                self.fits.append(kmtc_i_fit)
            elif i == 4:
                self.fits.append(kmtc_v_fit)
            else:
                if dataset.bandpass == 'I':
                    q_flux = (kmtc_i_fit.source_fluxes[1] /
                              kmtc_i_fit.source_fluxes[0])
                elif dataset.bandpass == 'V':
                    q_flux = (kmtc_v_fit.source_fluxes[1] /
                              kmtc_v_fit.source_fluxes[0])
                else:
                    raise Exception(
                        'Unknown bandpass: {0}. '.format(dataset.bandpass) +
                        'Fitting is only defined for I and V.')

                self.fits.append(
                    mm.FitData(model=self.model, dataset=dataset,
                               fix_q_flux=q_flux))
# RP - I haven't thought deeply about it, but this seems very complicated.

# Fit the fluxes
event = MyEvent(model=binary_source_model, datasets=datasets)
print(event.chi2)
print(event.source_fluxes)
print(event.blend_fluxes)
