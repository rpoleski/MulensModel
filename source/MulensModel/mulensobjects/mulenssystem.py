import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.constants import au, c, G

from MulensModel.mulensobjects.lens import Lens
from MulensModel.mulensobjects.source import Source
from MulensModel.model import Model


class MulensSystem(object):
    """
    A microlensing system consisting of a lens and a source.
    """

    def __init__(self, lens=None, source=None, mu_rel=None):
        self._lens = None
        self._source = None
        if lens is not None:
            if source is None:
                raise AttributeError(
                    'If lens is specified, source must also be specified.')
            else:
                self.lens = lens
                self.source = source
        if source is not None and lens is None:
            raise AttributeError(
                'If source is specified, lens must also be specified.')
        if mu_rel is not None:
            self.mu_rel = mu_rel
        else:
            self._mu_rel = None

    def __repr__(self):
        output_str = '------\n{0}\n{1}\ntheta_E = {2}\n'.format(
            self.lens, self.source, self.theta_E)
        if self._mu_rel is not None:
            mu_rel_str = 'mu_rel = {0}\n'.format(self._mu_rel)
            t_E_str = 't_E = {0}\n'.format(self.t_E)
            output_str += mu_rel_str + t_E_str
        return output_str + '------'

    @property
    def lens(self):
        """
        A :py:class:`~MulensModel.mulensobjects.lens.Lens` object.
        Physical properties of the lens. Note: lens mass must be in
        solMasses.
        """
        return self._lens

    @lens.setter
    def lens(self, value):
        if isinstance(value, Lens):
            self._lens = value
        else:
            raise TypeError("lens must be a Lens object")
        if (self.source is not None and self.source.distance is not None and
                self.lens.distance is not None and
                self.source.distance.value < self.lens.distance.value):
            msg = 'Source cannot be closer than lens: {:} {:}'
            raise ValueError(
                msg.format(self.source.distance, self.lens.distance))

    @property
    def source(self):
        """
        :py:class:`~MulensModel.mulensobjects.source.Source` object.
        Physical properties of the source.
        """
        return self._source

    @source.setter
    def source(self, value):
        if isinstance(value, Source):
            self._source = value
        else:
            raise TypeError(
                "source must be a MulensModel.mulensobjects.source.Source" +
                "object")
        if (self.lens is not None and self.source.distance is not None and
                self.lens.distance is not None and
                self.source.distance.value < self.lens.distance.value):
            msg = 'Source cannot be closer than lens: {:} {:}'
            raise ValueError(msg.format(
                self.source.distance, self.lens.distance))

    @property
    def mu_rel(self):
        """
        *astropy.Quantity*

        Relative proper motion between the source and lens
        stars. If set as a *float*, units are assumed to be mas/yr.
        """
        return self._mu_rel

    @mu_rel.setter
    def mu_rel(self, value):
        if isinstance(value, u.Quantity):
            self._mu_rel = value
        else:
            self._mu_rel = value * u.mas / u.yr

    @property
    def t_E(self):
        """
        *astropy.Quantity*

        The Einstein crossing time (in days). If set as a *float*,
        assumes units are in days.
        """
        try:
            t_E = self.theta_E/self.mu_rel
            return t_E.to(u.day)
        except Exception:
            return None

    @t_E.setter
    def t_E(self, t_E):
        if isinstance(t_E, u.Quantity):
            self.mu_rel = self.theta_E / t_E.to(u.year)
        else:
            self.mu_rel = self.theta_E / t_E * u.year

    @property
    def pi_rel(self):
        """
        *astropy.Quantity*, read-only

        The source-lens relative parallax in milliarcseconds.
        """
        return self.lens.pi_L.to(u.mas) - self.source.pi_S.to(u.mas)

    @property
    def pi_E(self):
        """
        *float*, read-only

        The Einstein ring radius. It's equal to pi_rel / theta_E.
        Dimensionless.
        """
        return (self.pi_rel / self.theta_E).decompose().value

    @property
    def theta_E(self):
        """
        *astropy.Quantity*, read-only

        The angular Einstein Radius in milliarcseconds.
        """
        kappa = (4. * G / (c**2 * au)).to(
            u.mas/u.Msun, equivalencies=u.dimensionless_angles())

        return np.sqrt(
            kappa * self.lens.total_mass.to(u.solMass) *
            self.pi_rel.to(u.mas))

    @property
    def r_E(self):
        """
        *astropy.Quantity*, read-only

        The physical size of the Einstein Radius in the Lens plane (in AU).
        """
        return (self.lens.distance * self.theta_E.to(
                '', equivalencies=u.dimensionless_angles())).to(u.au)

    @property
    def r_E_tilde(self):
        """
        *astropy.Quantity*, read-only

        The physical size of the Einstein Radius projected onto the
        Observer plane (in AU).
        """
        return self.r_E * self.source.distance / (
            self.source.distance - self.lens.distance)

    def plot_magnification(self, u_0=None, alpha=None, **kwargs):
        """
        Plot the magnification curve for the lens. u_0 must always be
        specified. If the lens has more than one body, alpha must also
        be specified.

        Parameters :
            u_0: *float*
                Impact parameter between the source and the lens (as a
                fraction of the Einstein ring)

            alpha: *astropy.Quantity*, *float*
                If *float* then degrees are assumed as a unit.
                See :py:obj:`MulensModel.modelparameters.ModelParameters.alpha`

            ``**kwargs``:
                See :py:func:`MulensModel.model.Model.plot_magnification()`
        """
        if u_0 is None:
            raise AttributeError('u_0 is required')
        else:
            parameters = {'t_0': 0., 'u_0': u_0}
            if self.t_E is not None:
                parameters['t_E'] = self.t_E
                xtitle = 'Time (days)'
            else:
                parameters['t_E'] = 1.
                xtitle = 'Time (tE)'

            if self.source.angular_radius is not None:
                parameters['rho'] = (self.source.angular_radius.to(u.mas) /
                                     self.theta_E.to(u.mas))

            if self.lens.n_masses > 1:
                parameters['q'] = self.lens.q
                parameters['s'] = self.lens.s
                if alpha is None:
                    raise AttributeError(
                        'alpha is required for 2-body lenses.')
                else:
                    parameters['alpha'] = alpha

            model = Model(parameters=parameters)
            model.plot_magnification(**kwargs)
            plt.xlabel(xtitle)

    def plot_caustics(self, n_points=5000, **kwargs):
        """
        Plot the caustics structure using `Pyplot scatter`_. See
        :py:func:`MulensModel.caustics.Caustics.plot()`

        Parameters :
            n_points: *int*
                Number of points be plotted.

            ``**kwargs``:
                Keyword arguments passed to `Pyplot scatter`

        .. _Pyplot scatter:
           https://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.scatter

        """
        self.lens.plot_caustics(n_points=n_points, **kwargs)
