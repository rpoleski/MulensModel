import numpy as np
import matplotlib.pyplot as pl

from astropy import units as u
from astropy.constants import au, c, G

from MulensModel.mulensobjects.lens import Lens
from MulensModel.mulensobjects.source import Source
from MulensModel.model import Model
from MulensModel.modelparameters import ModelParameters


class MulensSystem(object):
    """
    A microlensing system consisting of a lens and a source.
    """
    def __init__(self, lens=None, source=None, mu_rel=None):
        if lens is not None :
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
            output_str += mu_rel_str+t_E_str
        return output_str+'------'

    @property
    def lens(self):
        """
        Physical properties of the lens. A Lens object. Note: lens
        mass must be in solMasses.
        """
        return self._lens

    @lens.setter
    def lens(self, value):
        if isinstance(value, Lens):
            self._lens = value
        else:
            raise TypeError("lens must be a Lens object")

    @property
    def source(self):
        """
        Physical properties of the source. A Source object.
        """
        return self._source

    @source.setter
    def source(self, value):
        if isinstance(value, Source):
            self._source = value
        else:
            raise TypeError("source must be a Source object")

    @property
    def mu_rel(self):
        """
        Relative proper motion between the source and lens
        stars. (magnitude only)
        """
        return self._mu_rel

    @mu_rel.setter
    def mu_rel(self, value):
        if isinstance(value, u.Quantity):
            self._mu_rel = value
        else:
            self._mu_rel = value * u.mas / u.yr


    @property
    def pi_rel(self):
        """
        The source-lens relative parallax in milliarcseconds.
        """
        return self.lens.pi_L.to(u.mas) - self.source.pi_S.to(u.mas)

    @property
    def theta_E(self):
        """
        The angular Einstein Radius in milliarcseconds.
        """
        kappa = (4. * G / (c**2 * au)).to(u.mas/u.Msun, 
                        equivalencies=u.dimensionless_angles())

        return np.sqrt(
            kappa * self.lens.total_mass.to(u.solMass) 
            * self.pi_rel.to(u.mas))

    @property
    def r_E(self):
        """
        The physical size of the Einstein Radius in the Lens plane (in AU).
        """
        return (self.lens.distance * self.theta_E.to(
                '',equivalencies=u.dimensionless_angles())).to(u.au)

    @property
    def r_E_tilde(self):
        """
        The physical size of the Einstein Radius projected onto the
        Observer plane (in AU).
        """
        return self.r_E * self.source.distance / (self.source.distance - self.lens.distance)

    @property
    def t_E(self):
        """
        The Einstein crossing time (in days).
        """
        try:
            t_E = self.theta_E/self.mu_rel
            return t_E.to(u.day)
        except:
            return None

    @t_E.setter
    def t_E(self, t_E):
        if isinstance(t_E, u.Quantity):
            self.mu_rel = self.theta_E / t_E.to(u.year)
        else:
            self.mu_rel = self.theta_E / t_E * u.year

    def plot_magnification(self, u_0=None, alpha=None,**kwargs):
        """
        Plot the magnification curve for the lens. u_0 must always be
        specified. If the lens has more than one body, alpha must also
        be specified.
        """
        if u_0 is None:
            raise AttributeError('u_0 is required')
        else:
            parameters = ModelParameters(t_0=0., u_0=u_0)
            if self.t_E is not None:
                parameters.t_E = self.t_E
                xtitle = 'Time (days)'
            else:
                parameters.t_E = 1.
                xtitle = 'Time (tE)'

            if self.source.angular_size is not None:
                parameters.rho = (self.source.angular_size.to(u.mas) 
                                  / self.theta_E.to(u.mas))

            if self.lens.n_masses > 1:
                parameters.q = self.lens.q
                parameters.s = self.lens.s
                if alpha is None:
                    raise AttributeError('alpha is required for 2-body lenses.')
                else:
                    parameters.alpha = alpha

            model = Model(parameters=parameters)
            model.plot_magnification(**kwargs)
            pl.xlabel(xtitle)

    def plot_caustics(self, n_points=5000, **kwargs):
        """
        Plot the caustic structure. n_points specifies the number of
        points to generate in the caustic.
        """
        self.lens.plot_caustics(n_points=n_points, **kwargs)

