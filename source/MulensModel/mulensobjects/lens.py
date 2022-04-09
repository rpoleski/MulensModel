import numpy as np
from math import fsum
from astropy import units as u

from MulensModel.caustics import Caustics


class Lens(object):
    """
    A mass or system of masses.

    Standard parameter combinations for defining a Lens object:

    Point Lens:
        (:py:obj:`mass`, :py:obj:`distance`)

    2+ body lens:

        (:py:obj:`s`, :py:obj:`q`)

        (:py:obj:`epsilon`, :py:obj:`s`)

        (:py:obj:`total_mass`, :py:obj:`epsilon`, :py:obj:`distance`)
        - Missing s/Not Implemented?

        (:py:obj:`mass_1`, :py:obj:`mass_2`, :py:obj:`a_proj`,
        :py:obj:`distance`) - Not Implemented

    Note that s, q, and epsilon may be single values, lists, or numpy ndarrays.

    If units are not specified for a given mass, it is assumed the
    value given is in Solar Masses.

    If units are not specified for distance it is assumed the value is
    given in kpc.

    """

    def __init__(self, total_mass=None, mass=None, mass_1=None, mass_2=None,
                 a_proj=None, distance=None, q=None, s=None, epsilon=None):
        self._caustics = None

        # Set up system from s and q if defined
        if q is not None:
            self.q = q
            if s is not None:
                self.s = s
            else:
                msg = 'if q is defined, s must also be defined.'
                raise AttributeError(msg)

            # include system mass if appropriate
            if mass_1 is not None:
                self.total_mass = mass_1 * (1. + q)

        # Set total_mass and epsilon if defined
        if total_mass is not None:
            self.total_mass = total_mass
        if epsilon is not None:
            self._epsilon = np.array(epsilon)

        # Set total_mass, epsilon from mass_1 and mass_2
        if mass_2 is not None:
            self.total_mass = mass_1 + mass_2
            self._epsilon = np.array(
                [mass_1 / self.total_mass, mass_2 / self.total_mass])
        else:
            if mass_1 is not None and q is None:
                self._set_single_mass(mass_1)
            if mass is not None:
                self._set_single_mass(mass)

        # Set distance and projected separation
        if distance is not None:
            self.distance = distance
        else:
            self._distance = None
        if a_proj is not None:
            self._a_proj = a_proj

    def __repr__(self):
        """ Make a nice string representation of the object. """
        # Lens Distance
        try:
            dist_str = 'Lens Distance: {0}\n'.format(self._distance)
        except AttributeError:
            dist_str = ' '

        # Lens mass or q
        try:
            return(
                '{1}Lens Total Mass: {0}'.format(self._total_mass, dist_str))
        except AttributeError:
            return('{1}Lens mass ratio: {0}'.format(self.q, dist_str))
        else:
            return('Lens.py __repr__ error')

    @property
    def epsilon(self):
        """
        [*float, list, numpy.ndarray*]

        An array of mass fractions for each lens components:
        m_i/total_mass. Stored as a *numpy.ndarray*.

        """
        return self._epsilon

    @epsilon.setter
    def epsilon(self, new_epsilon):
        self._epsilon = np.array(new_epsilon)

    @property
    def q(self):
        """
        [*float, list, numpy.ndarray*]

        mass ratio for companions relative to the primary

        Array of mass ratios defined relative to the primary
        (m_i/m_1). Size is number of components -1. Set as a list,
        np.ndarray, or single value.

        Note: if :py:obj:`total_mass` is defined before q, it is
        assumed this is the mass of the primary. If you want this to
        actually be the total mass, define it after defining q.

        """
        if self._epsilon.size > 1:
            return self._epsilon[1:] / self._epsilon[0]
        else:
            raise AttributeError('Lens has only one body')

    @q.setter
    def q(self, new_q):
        # Update epsilon
        new_q = np.insert(new_q, 0, 1.)
        self._epsilon = new_q / fsum(new_q)

        try:
            if np.array(new_q).size == self._epsilon.size - 1:
                # Case 3: the entire lens is defined (new_q changes
                # the values of q)
                pass
            else:
                # Case 2: the primary is defined (new_q adds masses)
                if ((self._total_mass is not None) and
                        (self._last_mass_set != 'total_mass')):
                    self._total_mass = self._total_mass * fsum(new_q)
        except AttributeError:
            # Case 1: nothing is initialized (new_q directly sets epsilon)
            pass

    @property
    def s(self):
        """
        [*float, list, numpy.ndarray*]

        Separation between the components of the lens as a fraction of
        the Einstein ring. A np.ndarray or single value.

        Definitions for more than 2 lens bodies TBD
        """
        return self._s

    @s.setter
    def s(self, new_s):
        self._s = new_s

    @property
    def total_mass(self):
        """
        *astropy.Quantity*

        The total mass of the lens (sum of all components). An
        astropy.Quantity with mass units. If set as a *float*, units
        are assumed to be solMass.

        """
        return self._total_mass

    @total_mass.setter
    def total_mass(self, new_mass):
        if not isinstance(new_mass, u.Quantity):
            new_mass *= u.solMass
        elif new_mass.unit.physical_type == 'dimensionless':
            new_mass *= u.solMass
        elif new_mass.unit.physical_type != 'mass':
            msg = 'wrong physical_type of new total_mass: {:}'
            raise ValueError(msg.format(new_mass.unit.physical_type))

        self._total_mass = new_mass
        self._last_mass_set = 'total_mass'

    @property
    def mass(self):
        """
        *astropy.Quantity*

        The mass of a point lens --> total mass. An astropy.Quantity
        with mass units. May be set as a float (in which case solMass
        is assumed).
        """
        if self._epsilon.size > 1:
            raise TypeError(
                'mass can only be defined for a point lens. use total_mass' +
                'for multiple bodies')
        else:
            return self.total_mass

    @mass.setter
    def mass(self, new_mass):
        try:
            if self._epsilon.size > 1:
                raise TypeError('mass can only be defined for a point lens')
            else:
                self._set_single_mass(new_mass)
        except AttributeError:
            self._set_single_mass(new_mass)

    @property
    def mass_1(self):
        """
        *astropy.Quantity*

        The mass of the primary. Defined as total_mass *
        epsilon[0]. An *astropy.Quantity* with mass units. If set as a
        *float*, units are assumed to be solMass.

        """
        return self.total_mass * self._epsilon[0]

    @mass_1.setter
    def mass_1(self, new_mass):
        try:
            self._change_mass(new_mass, 0)
        except AttributeError:
            self._set_single_mass(new_mass)
        self._last_mass_set = 'mass_1'

    @property
    def mass_2(self):
        """
        *astropy.Quantity*

        The mass of the secondary. Defined as total_mass *
        epsilon[1]. An *astropy.Quantity* with mass units. If set as a
        *float*, units are assumed to be solMass.

        Note that if total_mass is defined before mass_2, and there is
        no epsilon corresponding to mass_2, mass_2 is added to the
        total_mass.

        """
        return self.total_mass * self._epsilon[1]

    @mass_2.setter
    def mass_2(self, new_mass, add=False):
        if self._epsilon.size > 1:
            self._change_mass(new_mass, 1)
        else:
            self._add_mass(new_mass, 1)
        self._last_mass_set = 'mass_1'

    @property
    def mass_3(self):
        """
        *astropy.Quantity*

        The mass of the tertiary. Defined as total_mass * epsilon[2].
        An *astropy.Quantity* with mass units. If set as a *float*,
        units are assumed to be solMass.

        Note that if total_mass is defined before mass_3, and there is
        no epsilon corresponding to mass_3, mass_3 is added to the total_mass.

        """
        return self.total_mass * self._epsilon[2]

    @mass_3.setter
    def mass_3(self, new_mass, add=False):
        if self._epsilon.size > 2:
            self._change_mass(new_mass, 2)
        else:
            self._add_mass(new_mass, 2)
        self._last_mass_set = 'mass_3'

    @property
    def n_masses(self):
        """
        *int*, read-only

        number of masses in the system.

        """
        try:
            return len(self._epsilon)
        except NameError:
            return 1
        else:
            return "lens.py: exception in Lens.n_masses"

    def _change_mass(self, new_mass, index):
        """
        Updates total_mass and epsilon array if the mass of one of the
        components is changed. e.g. mass_2 is changed from 1 MJup to 2
        MJup.

        """
        if not isinstance(new_mass, u.Quantity):
            new_mass *= u.solMass
        elif new_mass.unit.physical_type == 'dimensionless':
            new_mass *= u.solMass
        elif new_mass.unit.physical_type != 'mass':
            msg = 'wrong physical_type of new total_mass: {:}'
            raise ValueError(msg.format(new_mass.unit.physical_type))

        new_total_mass = self._total_mass(1. - self._epsilon[index]) + new_mass
        self._epsilon = self._total_mass * self._epsilon / new_total_mass
        self._total_mass = new_total_mass

    def _set_single_mass(self, new_mass):
        """
        Initializes total_mass and epsilon if only one mass componenet
        is defined (i.e. a point lens).
        """
        if isinstance(new_mass, u.Quantity):
            if new_mass.unit.physical_type == 'dimensionless':
                new_mass *= u.solMass
            elif new_mass.unit.physical_type != 'mass':
                msg = 'wrong physical_type of new total_mass: {:}'
                raise ValueError(msg.format(new_mass.unit.physical_type))
            self._total_mass = new_mass
        else:
            self._total_mass = new_mass * u.solMass

        self._epsilon = np.array([1.0])

    def _add_mass(self, new_mass, index):
        """
        Private function: Updates the total_mass and adds a component
        to the epsilon array if masses are added
        sequentially. e.g. the lens is defined by defining mass_1 and
        mass_2.
        """
        if not isinstance(new_mass, u.Quantity):
            new_mass *= u.solMass
        elif new_mass.unit.physical_type == 'dimensionless':
            new_mass *= u.solMass
        elif new_mass.unit.physical_type != 'mass':
            msg = 'wrong physical_type of new total_mass: {:}'
            raise ValueError(msg.format(new_mass.unit.physical_type))

        new_total_mass = self._total_mass + new_mass
        self._epsilon = self._total_mass * self._epsilon / new_total_mass
        self._epsilon = np.insert(
            self._epsilon, index, new_mass / new_total_mass)
        self._total_mass = new_total_mass

    @property
    def distance(self):
        """
        *astropy.Quantity*

        The distance to the lens.

        May be set as a *float*. If no unit is given, the value is
        assumed to be kpc.
        """
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        if not isinstance(new_distance, u.Quantity):
            self._distance = new_distance * 1000. * u.pc
        else:
            if new_distance.unit.physical_type != 'distance':
                TypeError('Wrong type of new_distance!')
            if (new_distance.unit == "pc") or (new_distance.unit == "kpc"):
                self._distance = new_distance
            else:
                raise u.UnitsError(
                    'Allowed units for Lens distance are "pc" or "kpc"')

    @property
    def pi_L(self):
        """
        *astropy.Quantity*

        The parallax to the lens in milliarcseconds. May be set as a
        *float*, in which case units are assumed to be
        milliarcseconds.

        """
        return self._distance.to(u.mas, equivalencies=u.parallax())

    @pi_L.setter
    def pi_L(self, new_value):
        if not isinstance(new_value, u.Quantity):
            new_value = new_value * u.mas
        self._distance = new_value.to(u.pc, equivalencies=u.parallax())

    @property
    def a_proj(self):
        """
        *astropy.Quantity*

        Projected separation between the components of the lens in
        AU. An *astropy.Quantity* with distance units. If set as *float*
        (without units), AU is assumed.

        """
        raise NotImplementedError('a_proj is not used, e.g. to set s')
        return self._a_proj

    @a_proj.setter
    def a_proj(self, new_a_proj):
        raise NotImplementedError('a_proj is not used, e.g. to set s')
        if not isinstance(new_distance, u.Quantity):
            new_a_proj = new_a_proj * u.au
        self._a_proj = new_a_proj

    @property
    def caustics(self):
        """
        A :py:class:`~MulensModel.caustics.Caustics` object, read-only
        """
        if self._caustics is None:
            if self.n_masses > 2:
                raise NotImplementedError(
                    'Caustics do not support more than 2 bodies')
            else:
                self._caustics = Caustics(q=self.q, s=self.s)
        return self._caustics

    def plot_caustics(self, n_points=5000, **kwargs):
        """
        A function to plot the x,y coordinates (scaled to the
        Einstein ring) of the caustics. `Pyplot scatter`_ is used for
        plotting. See :py:func:`MulensModel.caustics.Caustics.plot()`.

        Parameters :
            n_points: *int*
                Number of points be plotted.

            ``**kwargs``:
                Keyword arguments passed to `Pyplot scatter`

        .. _Pyplot scatter:
           https://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.scatter

        """
        self.caustics.plot(n_points=n_points, **kwargs)
