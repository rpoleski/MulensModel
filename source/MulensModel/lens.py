from astropy import units as u

message_max_masses = 'A maximum of 2 masses is supported.'

class Lens(object):
    def __init__(self, n_components=None,mass=None,mass_1=None,mass_2=None,
                 distance=None,q=None):
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
                self.n_components +=1
            self.mass_2 = mass_2
        if distance != None:
            self.distance = distance
        if q != None:
            while self.n_components < 2:
                self.n_components +=1
            self.q = q

    @property
    def total_mass(self):
        try:
            return self._total_mass
        except AttributeError:
            if self.n_components == 1:
                self._total_mass = self._mass_1
            elif self.n_components == 2:
                self._total_mass = self._mass_1+self._mass_2
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
        return self._mass_1

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
            return self._mass_2

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
                    self._mass_2 = new_mass
                    if self.n_components == 1:
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
        if self.n_components != 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            return self._mass_2/self._mass_1

    @q.setter
    def q(self,new_q):
        if self.n_components != 2:
            raise AttributeError('q is only defined if there are 2 bodies.')
        else:
            try:
                self._mass_2 = new_q*self._mass_1
            except:
                try:
                    self._mass_1 = self._total_mass/(1.+new_q) 
                    self._mass_2 = new_q*self._mass_1
                except AttributeError:
                    raise AttributeError('Either mass_1 or total_mass must be defined before q.')
