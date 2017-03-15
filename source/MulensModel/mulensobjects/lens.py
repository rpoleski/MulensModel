import numpy as np
from math import fsum
from astropy import units as u

class Lens(object):
    """
    A mass or system of masses.
    
    Attributes:
        total_mass: Total mass of the lens system
        epsilon: mass fraction for each component relative to the total mass
        q: mass ratio for companions relative to the primary
        s: separation between the companions and the primary as a
            fraction of the Einstein radius.
        mass: mass of the lens for a point mass
        mass_1: mass of the primary
        mass_2: mass of the secondary
        n_masses: total number of masses in the system
        distance: distance to the lens in pc.
        pi_L: parallax of the lens in mas.

    TO DO:
        - problem with tracking number of masses, esp when
          successively defining masses (see test_Lens.py)
        - implement triple+ systems
    """
    def __init__(self, total_mass=None, mass=None, mass_1=None, mass_2=None,
                 a_proj=None, distance=None, q=None, s=None, epsilon=None):
        """
        Define a Lens object. Standard parameter combinations:
        (s, q)
        (mass, distance)
        (mass_1, mass_2, a_proj, distance)
        (total_mass, epsilon, distance)
        (epsilon, s)
        Note that s, q, and epsilon may be lists or numpy arrays.

        If units are not specified for mass, it is assumed the value
        given is in Solar Masses. 

        If units are not specified for distance, and the value given
        is <50, it is assumed the value is given in kpc. Otherwise, pc
        are assumed.
        """
        self._caustic=None
        #Set up system from s and q if defined
        if q is not None:
            self.q = q
            if s is not None:
                self.s = s
            else:
                msg = 'if q is defined, s must also be defined.'
                raise AttributeError(msg)

            #include system mass if appropriate
            if mass_1 is not None:
                self.total_mass = mass_1 * (1. + q)

        #Set total_mass and epsilon if defined
        if total_mass is not None:
            self.total_mass = total_mass
        if epsilon is not None:
            self._epsilon = np.array(epsilon)

        #Set total_mass, epsilon from mass_1 and mass_2
        if mass_2 is not None:
            self.total_mass = mass_1 + mass_2
            self._epsilon = np.array(
                            [mass_1/self.total_mass, mass_2/self.total_mass])
        else:
            if mass_1 is not None and q is None:
                self._set_single_mass(mass_1)
            if mass is not None:
                self._set_single_mass(mass)

        #Set distance and projected separation
        if distance is not None:
            self.distance = distance
        if a_proj is not None:
            self._a_proj = a_proj

    def __repr__(self):
        """Make a nice string representation of the mass. NEEDS WORK."""
        dist_str = 'Lens Distance: {0}\n'.format(self._distance)
        try:
            return('{1}Lens Total Mass: {0}'.format(self._total_mass, dist_str))
        except NameError:
            return('{1}Lens components: {0}'.format(self._q, dist_str))
        else:
            return('Lens.py __repr__ error')

    @property
    def total_mass(self):
        """
        The total mass of the lens (sum of all components). An
        astropy.Quantity with mass units.
        """
        return self._total_mass

    @total_mass.setter
    def total_mass(self, new_mass):
        if not isinstance(new_mass, u.Quantity):
            new_mass = new_mass * u.solMass

        self._total_mass = new_mass

    @property
    def epsilon(self):
        """
        An array of mass ratios for each lens components:
        m_i/total_mass. A numpy array.
        """
        return self._epsilon

    @epsilon.setter
    def epsilon(self, new_epsilon):
        self._epsilon = np.array(new_epsilon)

    @property
    def n_masses(self):
        """
        number of masses in the system.
        """
        try:
            return len(self._epsilon)
        except NameError:
            return 1
        else:
            return "lens.py: exception in Lens.n_masses"

    @property
    def mass(self): 
        """
        The mass of a point lens --> total mass. An astropy.Quantity
        with mass units.
        """
        if self._epsilon.size > 1:
            raise TypeError(
                'mass can only be defined for a point lens. use total_mass for multiple bodies')
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
        The mass of the primary. Defined as total_mass *
        epsilon[0]. An astropy.Quantity with mass units.
        """
        return self.total_mass * self._epsilon[0]

    @mass_1.setter
    def mass_1(self, new_mass):
        try:
            self._change_mass(new_mass, 0)
        except AttributeError:
            self._set_single_mass(new_mass)

    @property
    def mass_2(self):
        """
        The mass of the secondary. Defined as total_mass *
        epsilon[1]. An astropy.Quantity with mass units.
        """
        return self.total_mass * self._epsilon[1]

    @mass_2.setter
    def mass_2(self, new_mass, add=False):
        """
        Note that if total_mass is defined before mass_2, and there is
        no epsilon corresponding to mass_2, this function will proceed
        to add mass_2 to the total_mass.
        """
        if self._epsilon.size > 1:
            self._change_mass(new_mass, 1)
        else:
            self._add_mass(new_mass, 1)

    @property
    def mass_3(self):
        """
        The mass of the tertiary. Defined as total_mass *
        epsilon[2]. An astropy.Quantity with mass units.
        """
        return self.total_mass * self._epsilon[2]

    @mass_3.setter
    def mass_3(self, new_mass, add=False):
        """
        Note that if total_mass is defined before mass_3, and there is
        no epsilon corresponding to mass_3, this function will proceed
        to add mass_3 to the total_mass.
        """
        if self._epsilon.size > 2:
            self._change_mass(new_mass, 2)
        else:
            self._add_mass(new_mass, 2)

    @property
    def distance(self):
        """
        The distance to the lens. An astropy Quantity.
        """
        return self._distance

    @distance.setter
    def distance(self, new_distance):
        """
        The distance should either be given in pc, or if no unit is
        given, the value is assumed to be kpc if it is <50 and in pc
        otherwise.
        """
        if not isinstance(new_distance, u.Quantity):
            if new_distance < 50:
                self._distance = new_distance * 1000. * u.pc
            else:
                self._distance = new_distance * u.pc
        else:
            if (new_distance.unit == "pc") or (new_distance.unit == "kpc"):
                self._distance = new_distance
            else:
                raise u.UnitsError(
                    'Allowed units for Lens distance are "pc" or "kpc"') 

    @property
    def pi_L(self):
        """
        The parallax to the lens in millarcseconds.
        """
        return self._distance.to(u.mas, equivalencies=u.parallax())

    @pi_L.setter
    def pi_L(self, new_value):
        if not isinstance(new_value, u.Quantity):
            new_value = new_value * u.mas
        self._distance = new_value.to(u.pc, equivalencies=u.parallax())

    @property
    def q(self):
        """
        Array of mass ratios defined relative to the primary (m_i/m_1). Size is
        number of components -1. A numpy.array or single value.
        """
        if self._epsilon.size > 1:
            return self._epsilon[1:] / self._epsilon[0]
        else:
            raise AttributeError('Lens has only one body')

    @q.setter
    def q(self, new_q):
        """
        Note: if total_mass is defined before q, it is assumed this is the
        mass of the primary. If you want this to actually be the total mass,
        define it after defining q.
        """
        #Update epsilon
        new_q = np.insert(new_q, 0, 1.)
        self._epsilon = new_q / fsum(new_q)

        #Update total_mass: DOES NOT LOOK LIKE IT WORKS RIGHT. Maybe
        #goes before update epsilon?
        try:
            if np.array(new_q).size == self._epsilon.size - 1:
        #Case 3: the entire lens is defined (new_q changes the values of q)
                pass
            else:
        #Case 2: the primary is defined (new_q adds masses)
                self._total_mass = self._total_mass * fsum(new_q)
        except AttributeError:
        #Case 1: nothing is initialized (new_q directly sets epsilon)
            pass

    @property
    def s(self):
        """
        Separation between the components of the lens as a fraction of
        the Einstein ring. A numpy.array or single value.
        """
        #Definitions for more than 2 lens bodies TBD
        return self._s

    @s.setter
    def s(self, new_s):
        self._s = new_s

    @property
    def a_proj(self):
        """
        Projected separation between the components of the lens in
        AU. An astropy.Quantity with distance units.
        """
        return self._a_proj

    @a_proj.setter
    def a_proj(self, new_a_proj):
        self._a_proj = new_a_proj

    def _change_mass(self, new_mass, index):
        """
        Private function: updates total_mass and epsilon array if the
        mass of one of the components is changed. e.g. mass_2 is
        changed from 1 MJup to 2 MJup.
        """
        if not isinstance(new_mass, u.Quantity):
            new_mass = new_mass * u.solMass

        new_total_mass = self._total_mass(1. - self._epsilon[index]) + new_mass
        self._epsilon = self._total_mass * self._epsilon / new_total_mass
        self._total_mass = new_total_mass

    def _set_single_mass(self, new_mass):
        """
        Private function: Initializes total_mass and epsilon if only
        one mass componenet is defined (i.e. a point lens).
        """
        if isinstance(new_mass, u.Quantity):
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
            new_mass = new_mass * u.solMass

        new_total_mass = self._total_mass + new_mass
        self._epsilon = self._total_mass * self._epsilon / new_total_mass
        self._epsilon = np.insert(self._epsilon, index, 
                                  new_mass / new_total_mass)
        self._total_mass = new_total_mass


    @property
    def caustic(self):
        """
        The x,y coordinates (scaled to the Einstein ring) of the
        caustic structure.
        """
        if self._caustic is None:
            self._caustic = self._calculate_caustics()
        return self._caustic

    def plot_caustics(self):
        """
        A function to plot the x,y coordinates (scaled to the
        Einstein ring) of the caustics.
        """
        raise NotImplementedError('please write Lens.plot_caustics')

    def _calculate_caustics(self):
        """
        Private function: calculate the x,y coordinates (scaled to the
        Einstein ring) of the caustics.
        """
        raise NotImplementedError('please write Lens._calculate_caustics')
