import warnings
import numpy as np
from math import fsum
from matplotlib import rcParams, gridspec
import matplotlib.pyplot as plt

from MulensModel.fitdata import FitData
from MulensModel.mulensdata import MulensData
from MulensModel.model import Model
from MulensModel.coordinates import Coordinates
from MulensModel.utils import PlotUtils


class Event(object):
    """
    Combines a microlensing model with data. Allows calculating chi^2 and
    making a number of plots.

    Arguments :
        :py:obj:`~datasets` :  :py:class:`~MulensModel.mulensdata.MulensData`
        or *list* of :py:class:`~MulensModel.mulensdata.MulensData` objects,
            Datasets that will be linked to the event. These datasets will
            be used for chi^2 calculation, plotting etc.

        :py:obj:`~model` : :py:class:`~MulensModel.model.Model`
            Microlensing model that will be linked to the event. In order to
            get chi^2 for different sets of model parameters you should
            keep a single :py:class:`~MulensModel.model.Model` instance and
            change parameters for this model (i.e., do not provide separate
            :py:class:`~MulensModel.model.Model` instances).

        :py:obj:`~coords` : *str*,
        :py:class:`~MulensModel.coordinates.Coordinates`, or astropy.SkyCoord_
            Coordinates of the event. If *str*, then needs format accepted by
            astropy.SkyCoord_ e.g., ``'18:00:00 -30:00:00'``.

        fix_blend_flux, fix_source_flux: *dict*
            Used to fix the source flux(es) or blend flux
            for a particular dataset. The dataset is
            the key, and the value to be fixed is the value. For example, to
            fix the blending of some dataset *my_data* to zero set
            *fix_blend_flux={my_data: 0.}*. See also
            :py:class:`~MulensModel.fitdata.FitData` .

        fix_source_flux_ratio: *dict*
            Used to fix the flux ratio for a given band or dataset. The keys
            should be either :py:class:`~MulensModel.mulensdata.MulensData`
            objects or *str*. If a
            :py:class:`~MulensModel.mulensdata.MulensData` object is specified,
            it will take precedence over a band.

        fit: DEPRECATED

        data_ref: *int* or :py:class:`~MulensModel.mulensdata.MulensData`
            Reference dataset. If *int* then gives index of reference dataset
            in :py:attr:`~datasets`. Default is the first dataset.

    The datasets can be in magnitude or flux spaces. When we calculate chi^2
    we do it in magnitude or flux space depending on value of
    :py:attr:`~MulensModel.mulensdata.MulensData.chi2_fmt` attribute.
    If dataset is in magnitude space and model results
    in negative flux, then we calculate chi^2 in flux space but only for the
    epochs with negative model flux.

    .. _astropy.SkyCoord:
      http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    """

    def __init__(
            self, datasets=None, model=None, coords=None, fix_blend_flux=None,
            fix_source_flux=None, fix_source_flux_ratio=None, data_ref=0):
        self._model = None
        self._coords = None

        # Initialize self._model (and check that model is defined).
        if isinstance(model, Model):
            self._model = model
        elif model is not None:
            raise TypeError('incorrect argument model of class Event()')

        # Initialize self._datasets (and check that datasets is defined).
        if isinstance(datasets, (list, tuple, MulensData)) or datasets is None:
            self._set_datasets(datasets)
        else:
            raise TypeError('incorrect argument datasets of class Event()')

        self._data_ref = self._set_data_ref(data_ref)

        # Set event coordinates
        if coords is not None:
            self._update_coords(coords=coords)
        elif self._model is not None:
            if self._model.coords is not None:
                self._update_coords(coords=self._model.coords)

        self.sum_function = 'numpy.sum'

        # Properties related to FitData
        self._fits = None  # New property
        self.chi2 = None
        if fix_blend_flux is None:
            self.fix_blend_flux = {}
        else:
            self.fix_blend_flux = fix_blend_flux

        if fix_source_flux is None:
            self.fix_source_flux = {}
        else:
            self.fix_source_flux = fix_source_flux

        if fix_source_flux_ratio is None:
            self.fix_source_flux_ratio = {}
        else:
            self.fix_source_flux_ratio = fix_source_flux_ratio

    def plot(self, t_range=None, residuals=True, show_errorbars=None,
             show_bad=None, legend=True, trajectory=None, title=None,
             subtract_2450000=True, subtract_2460000=False, data_ref=None):
        """
        Basic plotting. Default is to plot the light curve with the residuals.
        If the model has 2 or more lenses, also plot the source trajectory with
        caustics. For more detailed control over the plotting, see
        :py:func:`~plot_model()`, :py:func:`~plot_data()`, and
        :py:func:`~plot_trajectory`.

        Keywords:
            t_range: *list*, *tuple*
                Time range over which to show the light curve plot.

            residuals: *bool*
                Whether or not to plot the residuals with the light curve.

            show_errorbars: *bool*
                Whether or not to show the errorbars on the data. Default is
                *None* (see :py:func:`~plot_data()`).

            show_bad: *bool*
                Whether or not to show data points marked as "bad" (see
                :py:attr:`~MulensModel.mulensdata.MulensData.bad`). Default is
                *None* (see :py:func:`~plot_data()`).

            legend: *bool*
                Whether or not to show a legend for the datasets.

            trajectory: *bool*
                Whether or not to plot the source trajectory. Defaults to
                *False* for single lens events and *True* for binary lenses.

            title: *str*
                Title for the plot. Same title is used for trajectory plot, if
                applicable.

            subtract_2450000, subtract_2460000: *bool*
                see :py:func:`~plot_data()`.

            data_ref:  *int* or *MulensData*
                see :py:func:`~plot_data()`.
        """
        if trajectory is None:
            if self.model.n_lenses == 1:
                trajectory = False
            else:
                trajectory = True

        self._plot_lc_default(
            t_range=t_range, residuals=residuals,
            show_errorbars=show_errorbars,
            show_bad=show_bad, legend=legend, title=title,
            subtract_2450000=subtract_2450000,
            subtract_2460000=subtract_2460000, data_ref=data_ref)

        if trajectory:
            self._plot_trajectory_default(t_range=t_range, title=title)

    def _plot_lc_default(
             self, t_range=None, residuals=True, show_errorbars=None,
             show_bad=None, legend=True, title=None,
             subtract_2450000=False, subtract_2460000=False, data_ref=None):
        """
        Plot model and data
        """

        plt.figure()

        if residuals:
            gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
            if title is not None:
                plt.suptitle(title)

            ax11 = plt.subplot(gs[0])
        else:
            if title is not None:
                plt.title(title)

        self.plot_model(
            t_range=t_range, subtract_2450000=subtract_2450000,
            subtract_2460000=subtract_2460000, data_ref=data_ref,
            color='black')
        self.plot_data(
            show_errorbars=show_errorbars, show_bad=show_bad,
            subtract_2450000=subtract_2450000,
            subtract_2460000=subtract_2460000, data_ref=data_ref)
        plt.gca().minorticks_on()
        if legend:
            plt.legend()

        if residuals:
            plt.xlabel(None)
            # plt.gca().xaxis.set_ticklabels([])
            plt.subplot(gs[1], sharex=ax11)
            self.plot_residuals(
                subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000,
                show_errorbars=show_errorbars,
                show_bad=show_bad, data_ref=data_ref)
            plt.setp(ax11.get_xticklabels(), visible=False)
            plt.setp(ax11.get_xaxis().get_offset_text(), visible=False)

        if t_range is not None:
            if subtract_2450000:
                xlim = np.array(t_range) - 2450000.
            elif subtract_2460000:
                xlim = np.array(t_range) - 2460000.
            else:
                xlim = t_range

            plt.xlim(xlim)

        plt.tight_layout()

    def _plot_trajectory_default(self, t_range=None, title=None):
        """
        plot trajectory after plotting model and data
        """
        plt.figure()
        plt.gca().set_aspect('equal')
        if title is not None:
            plt.title(title)

        self.plot_trajectory(t_range=t_range, caustics=True)
        self.plot_source_for_datasets()

        if t_range is None:
            t_range = (self.model.parameters.t_0 +
                       self.model.parameters.t_E * np.array([-1., 1.]))

        lims = (np.array(t_range) -
                self.model.parameters.t_0) / self.model.parameters.t_E
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(lims)
        plt.ylim(lims)
        plt.gca().minorticks_on()
        plt.tight_layout()

    def plot_model(self, data_ref=None, **kwargs):
        """
        Plot the model light curve in magnitudes. See
        :py:func:`MulensModel.model.Model.plot_lc()` for details.

        Keywords :
            data_ref: *int* or *MulensData*
                If data_ref is not specified, uses :py:obj:`~data_ref`.

            ``**kwargs``:
                Keywords passed to
                :py:func:`MulensModel.model.Model.plot_lc()`.
                You can use them to set time range plotted, fluxes,
                limb-darkening etc.
        """
        if data_ref is None:
            data_ref = self.data_ref

        (f_source_0, f_blend_0) = self.get_flux_for_dataset(data_ref)
        self.model.plot_lc(
            source_flux=f_source_0, blend_flux=f_blend_0, **kwargs)

    def plot_data(
            self, phot_fmt='mag', data_ref=None, show_errorbars=None,
            show_bad=None,
            subtract_2450000=False, subtract_2460000=False, **kwargs):
        """
        Plot the data scaled to the model.

        Keywords (all optional):
            phot_fmt: *string* ('mag', 'flux')
                Whether to plot the data in magnitudes or in flux. Default
                is 'mag'.

            data_ref: *int* or *MulensData*
                If data_ref is not specified, uses :py:obj:`~data_ref`.

            show_errorbars: *boolean* or *None*
                Do you want errorbars to be shown for all datasets?
                Default is *None*, which means the option is taken from each
                dataset plotting properties (for which default is *True*).
                If *True*, then data are plotted using matplotlib.errorbar().
                If *False*, then data are plotted using matplotlib.scatter().

            show_bad: *boolean* or *None*
                Do you want data marked as bad to be shown?
                Default is *None*, which means the option is taken from each
                dataset plotting properties (for which default is *False*).
                If bad data are shown, then they are plotted with 'x' marker.

            subtract_2450000, subtract_2460000: *boolean*
                If True, subtracts 2450000 or 2460000 from the time
                axis to get more human-scale numbers. If using, make
                sure to also set the same settings for all other
                plotting calls (e.g. :py:func:`plot_lc()`).

            ``**kwargs``:
                Passed to matplotlib plotting functions. Contrary to
                previous behavior, ``**kwargs`` are no longer remembered.

        """
        self._set_default_colors()  # For each dataset
        if self.fits is None:
            self.get_chi2()

        if data_ref is None:
            data_ref = self.data_ref

        # JCY want to implement show_errobars, show_bad as list option, so it
        # can be different for different datasets. DO LATER.

        # Set plot limits
        t_min = 3000000.
        t_max = 0.
        subtract = PlotUtils.find_subtract(subtract_2450000, subtract_2460000)

        # Get fluxes for the reference dataset
        (f_source_0, f_blend_0) = self.get_flux_for_dataset(data_ref)
        for (i, data) in enumerate(self._datasets):
            # Scale the data flue
            (flux, err_flux) = self.fits[i].scale_fluxes(f_source_0, f_blend_0)
            (y_value, y_err) = PlotUtils.get_y_value_y_err(
                phot_fmt, flux, err_flux)

            data._plot_datapoints(
                (y_value, y_err), subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000,
                show_errorbars=show_errorbars, show_bad=show_bad, **kwargs)

            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))

        # Plot properties
        plt.ylabel('Magnitude')
        plt.xlabel(
            PlotUtils.find_subtract_xlabel(subtract_2450000, subtract_2460000))
        plt.xlim(t_min-subtract, t_max-subtract)

        (ymin, ymax) = plt.gca().get_ylim()
        if ymax > ymin:
            plt.gca().invert_yaxis()

    def plot_residuals(
            self, show_errorbars=None, data_ref=None, subtract_2450000=False,
            subtract_2460000=False, show_bad=None, **kwargs):
        """
        Plot the residuals (in magnitudes) to the model.

        Keywords:
            For explanation of keywords, see doctrings in
            :py:func:`plot_data()`. Note different order of keywords.

        """
        self._set_default_colors()

        if data_ref is None:
            data_ref = self.data_ref

        # Plot limit parameters
        t_min = 3000000.
        t_max = 0.
        subtract = PlotUtils.find_subtract(subtract_2450000, subtract_2460000)

        # Plot zeropoint line
        plt.plot([0., 3000000.], [0., 0.], color='black')

        # Plot residuals
        (f_source_0, f_blend_0) = self.get_flux_for_dataset(data_ref)
        for i, data in enumerate(self._datasets):
            # Evaluate whether or nor it is necessary to calculate the model
            # for bad datapoints.
            if show_bad:
                bad = True
            else:
                bad = False

            (residuals, errorbars) = self.fits[i].get_residuals(
                phot_fmt='scaled', source_flux=f_source_0,
                blend_flux=f_blend_0, bad=bad)
            y_value = residuals
            y_err = errorbars
            data._plot_datapoints(
                (y_value, y_err), subtract_2450000=subtract_2450000,
                subtract_2460000=subtract_2460000,
                show_errorbars=show_errorbars, show_bad=show_bad, **kwargs)

            t_min = min(t_min, np.min(data.time))
            t_max = max(t_max, np.max(data.time))

        # Plot properties
        y_lim = np.max([np.abs(y_lim) for y_lim in plt.gca().get_ylim()])
        if y_lim > 1.:
            y_lim = 0.5

        plt.ylim(y_lim, -y_lim)
        plt.xlim(t_min-subtract, t_max-subtract)
        plt.ylabel('Residuals')
        plt.xlabel(
            PlotUtils.find_subtract_xlabel(subtract_2450000, subtract_2460000))

    def plot_trajectory(self, **kwargs):
        """
        Plot the trajectory of the source. See
        :py:func:`MulensModel.model.Model.plot_trajectory()` for details.
        """
        if 'show_data' in kwargs:
            raise AttributeError('Parameter show_data is deprecated. Use '
                                 'plot_source_for_datasets() instead.')

        self.model.plot_trajectory(**kwargs)

    def plot_source_for_datasets(self, **kwargs):
        """
        Plot source positions for all linked datasets.
        See :py:func:`MulensModel.model.Model.plot_source()` for
        details.

        Note: plots all points in datasets (including ones flagged as bad)
        using the same marker.
        """
        self._set_default_colors()

        for dataset in self.datasets:
            # RP: Call to MulensData private function should be replaced by
            # public functions.
            properties = dataset._set_plot_properties()
            self.model.plot_source(
                times=dataset.time, color=properties['color'], **kwargs)

    def _set_default_colors(self):
        """
        If the user has not specified a color for a dataset, assign
        one.
        """
        colors = [cycle['color'] for cycle in rcParams['axes.prop_cycle']]

        # Below we change the order of colors to most distinct first.
        used_colors = []
        for data in self._datasets:
            if 'color' in data.plot_properties.keys():
                used_colors.append(data.plot_properties['color'])

        if len(used_colors) == len(self._datasets):
            return

        if len(used_colors) == 0:
            differences = None
        else:
            diffs = np.array(
                [np.min(
                    PlotUtils.get_color_differences(used_colors, c))
                 for c in colors])
            indexes = np.argsort(diffs)[::-1]
            colors = [colors[i] for i in indexes]
            differences = diffs[indexes]

        # Assign colors when needed.
        color_index = 0
        for data in self._datasets:
            if 'color' not in data.plot_properties.keys():
                if differences is not None:
                    if differences[color_index] < 0.35:
                        msg = ('The color assign to one of the datasets in ' +
                               'automated way (' + colors[color_index] +
                               ') is very similar to already used color')
                        warnings.warn(msg, UserWarning)

                data.plot_properties['color'] = colors[color_index]
                color_index += 1
                if color_index == len(colors):
                    color_index = 0
                    msg = ('Too many datasets without colors assigned - ' +
                           'same color will be used for different datasets')
                    warnings.warn(msg, UserWarning)

    def get_flux_for_dataset(self, dataset):
        """
        Get the source and blend flux for a given dataset.

        Parameters :
            dataset: :py:class:`~MulensModel.mulensdata.MulensData` or *int*
            If *int* should be the index (starting at 0) of the appropriate
            dataset in the :py:obj:`~datasets` list.

        Returns :
            source_flux: *np.ndarray*
                flux of sources. see
                :py:obj:`~MulensModel.fitdata.FitData.source_fluxes`
            blend_flux: *float*
                blending flux. see
                :py:obj:`~MulensModel.fitdata.FitData.blend_flux`

        NOTE: This function does not recalculate fits or fluxes. If the data
        haven't yet been fit to the model (i.e. self.fits = None),
        it will run :py:func:`~fit_fluxes()`. Otherwise, it just accesses the
        existing values. So if you change something in :py:obj:`~model` or
        some fit parameter (e.g., :py:obj:`~fix_blend_flux`), be sure to run
        :py:func:`~fit_fluxes()` first.

        """
        if self.fits is None:
            self.fit_fluxes()

        if isinstance(dataset, MulensData):
            i = self.datasets.index(dataset)
        else:
            i = dataset

        source_flux = self.fits[i].source_fluxes
        blend_flux = self.fits[i].blend_flux

        return (source_flux, blend_flux)

    def get_ref_fluxes(self, data_ref=None, fit_blending=None):
        """
        Get source and blending fluxes for the reference dataset. See
        :py:func:`~get_flux_for_dataset()`. If the reference dataset is not
        set, uses the first dataset as default. See :py:obj:`~data_ref`.
        """
        if data_ref is not None:
            warnings.warn(
                'data_ref will be deprecated. It is redundant for getting ' +
                'the flux of the reference dataset. For the flux of an ' +
                'arbitrary dataset, use get_flux_for_dataset')

        if fit_blending is not None:
            self._apply_fit_blending(fit_blending)

        return self.get_flux_for_dataset(self.data_ref)

    def get_chi2(self, fit_blending=None):
        """
        Calculates chi^2 of current model by fitting for source and
        blending fluxes.

        Parameters :
            fit_blending: DEPRECATED. use :py:attr:`~fix_blend_flux` instead.

        Returns :
            chi2: *float*
                Chi^2 value

        """
        if fit_blending is not None:
            self._apply_fit_blending(fit_blending)

        self.fit_fluxes()
        chi2 = []
        for (i, dataset) in enumerate(self.datasets):
            # Calculate chi2 for the dataset excluding bad data
            chi2.append(self.fits[i].chi2)

        self.chi2 = self._sum(chi2)

        return self.chi2

    def get_chi2_for_dataset(self, index_dataset, fit_blending=None):
        """
        Calculates chi^2 for a single dataset

        Parameters :
            index_dataset: *int*
                index that specifies for which dataset the chi^2 is requested

            fit_blending: DEPRECATED. use :py:attr:`~fix_blending` instead.

        Returns :
            chi2: *float*
                chi2 for dataset[index_dataset].

        """
        if fit_blending is not None:
            self._apply_fit_blending(fit_blending)

        self.fit_fluxes()

        return self.fits[index_dataset].chi2

    def get_chi2_per_point(self, fit_blending=None, bad=False):
        """
        Calculates chi^2 for each data point of the current model by
        fitting for source and blending fluxes.

        Parameters :
            fit_blending: DEPRECATED. use :py:attr:`~fix_blending` instead.

            bad: *bool*
                Should chi2 be also caclulated for points marked as bad in
                MulensData? If `False` (default), then bad epochs have chi2 of
                `np.nan`.

        Returns :
            chi2: *list* of *np.ndarray*
                Chi^2 contribution from each data point,
                e.g. ``chi2[data_num][k]`` returns the chi2 contribution
                from the *k*-th point of dataset *data_num*.

        Example :
            Assuming ``event`` is instance of Event class to get chi2
            for 10-th point point of 0-th dataset.

            .. code-block:: python

               chi2 = event.get_chi2_per_point()
               print(chi2[0][10])

        """
        if fit_blending is not None:
            self._apply_fit_blending(fit_blending)

        self.fit_fluxes(bad=bad)

        # Calculate chi^2 given the fit
        chi2_per_point = []
        for (i, dataset) in enumerate(self.datasets):
            chi2 = self.fits[i].chi2_per_point
            if not bad:
                chi2[dataset.bad] = np.nan
            chi2_per_point.append(chi2)

        return chi2_per_point

    def get_chi2_gradient(self, parameters, fit_blending=None):
        """
        Fit for fluxes and calculate chi^2 gradient (also called Jacobian),
        i.e., :math:`d chi^2/d parameter`.

        Parameters :
            parameters: *str* or *list*, required
                Parameters with respect to which gradient is calculated.
                Currently accepted parameters are: ``t_0``, ``u_0``, ``t_eff``,
                ``t_E``, ``pi_E_N``, and ``pi_E_E``. The parameters for
                which you request gradient must be defined in py:attr:`~model`.

            fit_blending: DEPRECATED. use :py:attr:`~fix_blending` instead.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient

        """
        if fit_blending is not None:
            self._apply_fit_blending(fit_blending)

        self.fit_fluxes()
        self.calculate_chi2_gradient(parameters)
        return self.chi2_gradient

    def calculate_chi2_gradient(self, parameters):
        """
        Calculate chi^2 gradient (also called Jacobian), i.e.,
        :math:`d chi^2/d parameter`.

        Parameters :
            parameters: *str* or *list*, required
                Parameters with respect to which gradient is calculated.
                Currently accepted parameters are: ``t_0``, ``u_0``, ``t_eff``,
                ``t_E``, ``pi_E_N``, and ``pi_E_E``. The parameters for
                which you request gradient must be defined in py:attr:`~model`.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient

        NOTE: Because this is not a 'get' function, it ASSUMES you have ALREADY
        fit for the fluxes, e.g. by calling get_chi2().
        """
        gradient = {param: 0 for param in parameters}
        for i, dataset in enumerate(self.datasets):
            data_gradient = self.fits[i].calculate_chi2_gradient(parameters)
            for j, p in enumerate(parameters):
                gradient[p] += data_gradient[j]

        if len(parameters) == 1:
            out = gradient[parameters[0]]
        else:
            out = np.array([gradient[p] for p in parameters])

        self._chi2_gradient = out

        return self._chi2_gradient

    def fit_fluxes(self, bad=False):
        """
        Fit for the optimal fluxes for each dataset (and its chi2)
        """

        self._fits = []
        for dataset in self.datasets:
            if dataset in self.fix_blend_flux.keys():
                fix_blend_flux = self.fix_blend_flux[dataset]
            else:
                fix_blend_flux = False

            if dataset in self.fix_source_flux.keys():
                fix_source_flux = self.fix_source_flux[dataset]
            else:
                fix_source_flux = False

            # JCY - This needs a unit test.
            if dataset in self.fix_source_flux_ratio.keys():
                fix_source_flux_ratio = self.fix_source_flux_ratio[dataset]
            else:
                if dataset.bandpass in self.fix_source_flux_ratio.keys():
                    fix_source_flux_ratio = self.fix_source_flux_ratio[
                        dataset.bandpass]
                else:
                    fix_source_flux_ratio = False

            fit = FitData(
                model=self.model, dataset=dataset,
                fix_blend_flux=fix_blend_flux, fix_source_flux=fix_source_flux,
                fix_source_flux_ratio=fix_source_flux_ratio)
            fit.update(bad=bad)  # Fit the fluxes and calculate chi2.
            self.fits.append(fit)

    def _sum(self, data):
        """calculate sum of the data"""
        if self.sum_function == 'numpy.sum':
            return np.sum(data)
        elif self.sum_function == 'math.fsum':
            return fsum(data)
        else:
            raise ValueError(
                'Event.sum_function unrecognized: ' + self.sum_function)

    @property
    def coords(self):
        """
        see :py:class:`~MulensModel.coordinates.Coordinates`
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        self._update_coords(coords=new_value)

    def _update_coords(self, coords=None):
        """Set the coordinates as a SkyCoord object"""
        self._coords = Coordinates(coords)

        if self._model is not None:
            self._model.coords = self._coords

        # We run the command below with try, because _update_coords() is called
        # by _set_datasets before self._datasets is set.
        try:
            for dataset in self._datasets:
                dataset.coords = self._coords
        except Exception:
            pass

    @property
    def model(self):
        """an instance of :py:class:`~MulensModel.model.Model`"""
        return self._model

    @model.setter
    def model(self, new_value):
        if not isinstance(new_value, Model):
            raise TypeError((
                'wrong type of Event.model: {:} instead of ' +
                'MulensModel.Model()').format(type(new_value)))
        self._model = new_value

        if new_value.coords is not None:
            self._update_coords(coords=new_value.coords)

        self._fits = None  # reset the fits if the model changed.

    @property
    def datasets(self):
        """
        a *list* of :py:class:`~MulensModel.mulensdata.MulensData`
        instances.
        """
        return self._datasets

    @datasets.setter
    def datasets(self, new_value):
        self._set_datasets(new_value)

    def _set_datasets(self, new_value):
        """
        sets the value of self._datasets
        can be called by __init__ or @datasets.setter
        passes datasets to property self._model
        """
        if isinstance(new_value, list):
            for dataset in new_value:
                if dataset.coords is not None:
                    self._update_coords(coords=dataset.coords)

        if isinstance(new_value, MulensData):
            if new_value.coords is not None:
                self._update_coords(coords=new_value.coords)

            new_value = [new_value]

        if new_value is None:
            self._datasets = None
            return

        self._datasets = new_value
        self._fits = None  # reset the fits if the data changed

    @property
    def data_ref(self):
        """
        Reference data set for scaling the model fluxes to (for
        plotting). May be set as a
        :py:class:`~MulensModel.mulensdata.MulensData` object or an
        index (*int*). Default is the first data set.

        Returns :
            index (*int*) of the relevant dataset.
        """
        if self._data_ref is None:
            return 0
        else:
            return self._data_ref

    @data_ref.setter
    def data_ref(self, new_value):
        self._set_data_ref(new_value)

    def _set_data_ref(self, new_value):
        """
        Set reference dataset. Not covered by unit tests.
        """
        if isinstance(new_value, MulensData):
            index = self.datasets.index(new_value)
            try:
                ind_2 = self.datasets.index(new_value, index+1)
            except ValueError:
                pass
            else:
                raise ValueError(
                    'Dataset is included in Event.datasets more than once.')

            self._data_ref = index
        elif isinstance(new_value, (int, np.int_)):
            self._data_ref = new_value
        else:
            raise TypeError(
                'data_ref must be set using either *int* or *MulensData*: ' +
                '{0}'.format(type(new_value)))

    @property
    def chi2(self):
        """
        *float*

        Chi^2 value. Note this is a static property. It is only updated when
        :py:func:`~fit_fluxes()` or :py:func:`~get_chi2()` is run. So, if you
        change one of the settings be sure to run one of those functions to
        update the chi2.
        """
        return self._chi2

    @chi2.setter
    def chi2(self, new_value):
        self._chi2 = new_value

    @property
    def chi2_gradient(self):
        """
        Return previously calculated chi^2 gradient (also called Jacobian),
        i.e., :math:`d chi^2/d parameter`. See :py:func:`~get_chi2_gradient()`
        and :py:func:`~calculate_chi2_gradient()`.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient. Will return None if the chi2 gradient was not
                previously calculated using one of the functions mentioned
                above.

        """
        try:
            return self._chi2_gradient
        except AttributeError:
            return None

    @property
    def fits(self):
        """
        *list* of :py:class:`~MulensModel.fitdata.FitData` objects

        There is one :py:class:`~MulensModel.fitdata.FitData` object for
        each dataset containing the information for fitting the model to
        that dataset, e.g. fitted fluxes, chi2 (for that dataset).
        """
        return self._fits

    @property
    def fluxes(self):
        """
        *list*

        An array giving the fitted source and blend flux(es) for each dataset.
        """
        fluxes = []
        if self.fits is None:
            raise AttributeError(
                'Fluxes not calculated. Run fit_fluxes() or get_chi2() first.')

        for fit in self.fits:
            fluxes.append([fit.source_fluxes, fit.blend_flux])

        return fluxes

    @property
    def source_fluxes(self):
        """
        *list*

        An array giving the fitted source flux(es) for each dataset.
        """
        fluxes = []
        if self.fits is None:
            raise AttributeError(
                'Fluxes not calculated. Run fit_fluxes() or get_chi2() first.')

        for fit in self.fits:
            fluxes.append(fit.source_fluxes)

        return fluxes

    @property
    def blend_fluxes(self):
        """
        *list*

        An array giving the fitted blend flux for each dataset.
        """
        fluxes = []
        if self.fits is None:
            raise AttributeError(
                'Fluxes not calculated. Run fit_fluxes() or get_chi2() first.')

        for fit in self.fits:
            fluxes.append(fit.blend_flux)

        return fluxes

    @property
    def sum_function(self):
        """
        *str*

        Function used for adding chi^2 contributions. Can be either
        'numpy.sum' (default value) or 'math.fsum'.
        The latter is slightly slower and more accurate,
        which may be important for large datasets.
        """
        return self._sum_function

    @sum_function.setter
    def sum_function(self, new_value):
        self._sum_function = new_value

    # ----Stuff that Doesn't Work (or is Deprecated)---- #
    def reset_best_chi2(self):
        """
        DEPRECATED

        Reset :py:attr:`~best_chi2` attribute and its parameters
        (:py:attr:`~best_chi2_parameters`).
        """
        raise AttributeError(
            'reset_best_chi2 (and best_chi2) has been deprecated.')

    @property
    def best_chi2(self):
        """
        DEPRECATED

        *float*

        The smallest value returned by :py:func:`get_chi2()`.
        """
        raise AttributeError('best_chi2 has been deprecated.')

    @property
    def best_chi2_parameters(self):
        """
        DEPRECATED

        *dict*

        Parameters that gave the smallest chi2.
        """
        raise AttributeError(
            'best_chi2_parameters (and best_chi2) has been deprecated.')

    def _apply_fit_blending(self, fit_blending):
        warnings.warn(
            'fit_blending option will be deprecated in future.' +
            'To fix the blending, set Event.fix_blend_flux instead.',
            FutureWarning)
        self._fits = None
        if fit_blending is True:
            self.fix_blend_flux = {}
        else:
            for dataset in self.datasets:
                self.fix_blend_flux[dataset] = 0.
