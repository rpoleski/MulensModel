from astropy import units as u

message_max_masses = 'A maximum of 2 masses is supported.'

class Lens(object):
    def __init__(self, n_components=None,mass=None,mass_1=None,mass_2=None,
                 a_proj=None,distance=None,q=None,s=None):
        if n_components > 2:
            raise ValueError(message_max_masses)
        else:
            self.n_components = n_components
        if mass != None:
            if n_components == None:
                self.n_components = 1
            self.mass = mass
        if mass_1 != None:
            if n_components == None:
                self.n_components = 1
            self.mass_1 = mass_1
        if mass_2 != None:
            while self.n_components < 2:
                try:
                    self.n_components +=1
                except AttributeError:
                    self.n_components = 2
            self.mass_2 = mass_2
        if distance != None:
            self.distance = distance
        if q != None:
            self.n_components = 2
            self.q = q
        if s != None:
            self.s = s
        if a_proj != None:
            self.a_proj = a_proj

    def __repr__(self):
        pass

    @property
    def total_mass(self):
        try:
            return self._total_mass
        except AttributeError:
            if self.n_components == 1:
                self._total_mass = self._mass_1
            elif self.n_components == 2:
                self._total_mass = self._mass_1*(1.+self._q)
            else:
                raise ValueError(message_max_masses)
            return self._total_mass

    @total_mass.setter
    def total_mass(self,new_mass):
        if type(new_mass) != u.Quantity:
            raise TypeError('mass_1 must have astropy units, i.e. it must be an astropy.units Quantity.')
        else:
            if (new_mass.unit == "solMass" or new_mass.unit == "jupiterMass"
                or new_mass.unit == "earthMass"):
                self._total_mass = new_mass
            else:
                raise u.UnitsError('Allowed units are "solMass", "jupiterMass", or "earthMass"')


    @property
    def mass(self):
        if self.n_components == 1:
            return self._mass_1
        else:
            return self.total_mass

    @mass.setter
    def mass(self,new_mass):
        if type(new_mass) != u.Quantity:
            raise TypeError('mass_1 must have astropy units, i.e. it must be an astropy.units Quantity.')
        else:
            if (new_mass.unit == "solMass" or new_mass.unit == "jupiterMass"
                or new_mass.unit == "earthMass"):
                self._mass_1 = new_mass
            else:
                raise u.UnitsError('Allowed units are "solMass", "jupiterMass", or "earthMass"')

    @property
    def mass_1(self):
        try:
            return self._mass_1
        except AttributeError:
            return self._total_mass/(1.+self._q)

    @mass_1.setter
    def mass_1(self,new_mass):
        if type(new_mass) != u.Quantity:
            raise TypeError('mass_1 must have astropy units, i.e. it must be an astropy.units Quantity.')
        else:
            if (new_mass.unit == "solMass" or new_mass.unit == "jupiterMass"
                or new_mass.unit == "earthMass"):
                self._mass_1 = new_mass
            else:
                raise u.UnitsError('Allowed units are "solMass", "jupiterMass", or "earthMass"')

    @property
    def mass_2(self):
        if self.n_components == 1:
            raise ValueError('Only one mass is defined')
        else:
            try:
                return self._q*self._mass_1
            except AttributeError:
                return self._total_mass/((1./self._q)+1.)

    @mass_2.setter
    def mass_2(self,new_mass,add=False):
        if self.n_components == 1 and add == False:
            raise ValueError('n_components=1; either increase the number of components or use add=True.')
        else:
            if type(new_mass) != u.Quantity:
                raise TypeError('mass_1 must have astropy units, i.e. it must be an astropy.units Quantity.')
            else:
                if (new_mass.unit == "solMass" 
                    or new_mass.unit == "jupiterMass"
                    or new_mass.unit == "earthMass"):
                    try: 
                        self._q = (new_mass/self._mass_1).decompose()
                    except:
                        try:
                            self._q = (new_mass
                                       /(self._total_mass-new_mass)).decompose()
                            self._mass_1 = new_mass/self._q
                        except AttributeError:
                            raise AttributeError('Either mass_1 or total_mass must be defined before mass_2.')
                    if self.n_components < 2:
                        self.n_components += 1
                else:
                    raise u.UnitsError('Allowed units are "solMass", "jupiterMass", or "earthMass"')


    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self,new_distance):
        if type(new_distance) != u.Quantity:
            raise TypeError('distance must have astropy units, i.e. it must be an astropy.units Quantity.')
        else:
            if (new_distance.unit == "pc" or new_distance.unit == "AU" 
                or new_distance.unit == "lyr"):
                self._distance = new_distance
            else:
                raise u.UnitsError('Allowed units are "pc", "AU", or "lyr"') 

    @property
    def q(self):
        if self.n_components < 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            return self._q

    @q.setter
    def q(self,new_q):
        if self.n_components < 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            self._q = new_q

    @property
    def s(self):
        if self.n_components < 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            return self._s

    @s.setter
    def s(self,new_s):
        if self.n_components < 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            self._s = new_s

    @property
    def a_proj(self):
        pass

    @a_proj.setter
    def a_proj(self,new_a_proj):
        pass
