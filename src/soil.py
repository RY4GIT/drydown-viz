# %%
"""
This code is based on the following source: 

natasha, Kelly Caylor, Noah Spahn, & Bryn Morgan. (2021). ecohydro/maize-Toff: First release of maize-Toff (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5204423

"""

# %%
from numpy import log as ln
from numpy import exp

#%% Set parameters related to soils and siginificant digits
rho = 1000.0    # density of water in kg/m^3
g = 9.8         # acceleration of gravity in m/s^2
field_capacity = -33 / 1000 # Field capacity in MPa.
PRECISION = 2   # Number of decimal places of precision in calculations (default is 2)

# DATA FROM CLAPP AND HORBERGER (C&H) 1978, Table 2:
soils = {
    'sand':{
        'b': 4.05,
        'Psi_S_cm': 12.1,   # saturated water tension, cm
        'Psi_l_cm': 4.66,   # Drainage water tension, cm
        'n': 0.395,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 1.056,        # saturated hydraulic conductivity, cm/min
        'S': 1.52           # sorptivity, cm/min^1/2    
    },
    'loamy sand':{
        'b': 4.38,
        'Psi_S_cm': 9.0,    # saturated water tension, cm
        'Psi_l_cm': 2.38,   # Drainage water tension, cm
        'n': 0.410,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.938,        # saturated hydraulic conductivity, cm/min
        'S': 1.04           # sorptivity, cm/min^1/2  
    },
    'sandy loam':{
        'b': 4.90,
        'Psi_S_cm': 21.8,   # saturated water tension, cm
        'Psi_l_cm': 9.52,   # Drainage water tension, cm
        'n': 0.435,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.208,        # saturated hydraulic conductivity, cm/min
        'S': 1.03           # sorptivity, cm/min^1/2  
    },
    'silt loam':{
        'b': 5.30,
        'Psi_S_cm': 78.6,   # saturated water tension, cm
        'Psi_l_cm': 75.3,   # Drainage water tension, cm
        'n': 0.485,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0432,       # saturated hydraulic conductivity, cm/min
        'S': 1.26           # sorptivity, cm/min^1/2  
    },
    'loam':{
        'b': 5.39,
        'Psi_S_cm': 47.8,   # saturated water tension, cm
        'Psi_l_cm': 20.0,   # Drainage water tension, cm
        'n': 0.451,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0417,       # saturated hydraulic conductivity, cm/min
        'S': 0.693          # sorptivity, cm/min^1/2  
    },
    'sandy clay loam':{
        'b': 7.12,
        'Psi_S_cm': 29.9,   # saturated water tension, cm
        'Psi_l_cm': 11.7,   # Drainage water tension, cm
        'n': 0.420,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0378,       # saturated hydraulic conductivity, cm/min
        'S': 0.488          # sorptivity, cm/min^1/2  
    },
    'silty clay loam':{
        'b': 7.75,
        'Psi_S_cm': 35.6,   # saturated water tension, cm
        'Psi_l_cm': 19.7,   # Drainage water tension, cm
        'n': 0.477,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0102,       # saturated hydraulic conductivity, cm/min
        'S': 0.310          # sorptivity, cm/min^1/2  
    },
    'clay loam':{
        'b': 8.52,
        'Psi_S_cm': 63.0,   # saturated water tension, cm
        'Psi_l_cm': 48.1,   # Drainage water tension, cm
        'n': 0.476,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0147,       # saturated hydraulic conductivity, cm/min
        'S': 0.537          # sorptivity, cm/min^1/2  
    },
    'sandy clay':{
        'b': 10.4,
        'Psi_S_cm': 15.3,   # saturated water tension, cm
        'Psi_l_cm': 8.18,   # Drainage water tension, cm
        'n': 0.426,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0130,       # saturated hydraulic conductivity, cm/min
        'S': 0.223          # sorptivity, cm/min^1/2  
    },
    'silty clay':{
        'b': 10.4,
        'Psi_S_cm': 49.0,   # saturated water tension, cm
        'Psi_l_cm': 23.0,   # Drainage water tension, cm
        'n': 0.492,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0062,       # saturated hydraulic conductivity, cm/min
        'S': 0.242          # sorptivity, cm/min^1/2  
    },
    'clay':{
        'b': 11.4,
        'Psi_S_cm': 40.5,   # saturated water tension, cm
        'Psi_l_cm': 24.3,   # Drainage water tension, cm
        'n': 0.482,         # porosity, cm^3/cm^3 (is Psi_S) in C&H,
        'Ks': 0.0077,       # saturated hydraulic conductivity, cm/min
        'S': 0.268          # sorptivity, cm/min^1/2  
    }
}

#%% Soil CLass Definition

class Soil():
    """ Defines a soil object based on either passed parameters
    or a soil texture class corresponding to the textures defined in 
    Clapp & Hornberger (C&H), 1978, Table 2.

    Usage: Soil(texture, params)
        Notes: 
            If texture is not provided, params must be included.
            If texture is provided, params is ignored.
            Capitilization in texture classes is ignored.
        
        texture = texture name from Clapp & Hornberger, 1978, Table 2.
            Valid options are:
            ["Sand", "Loamy Sand", "Sandy Loam", "Silt Loam", "Loam",
             "Sandy Clay Loam", Silty Clay Loam", "Clay Loam", "Sandy Clay",
             "Silty Clay", "Clay"]
        
        params = dictionary containing values for the soil parameters:
            params = {
                'b': 11.4,
                'Psi_S': 40.5,  # saturated water tension, cm
                'Psi_l': 24.3,  # Drainage water tension, cm
                'n': 0.482,     # porosity, cm^3/cm^3 (is Psi_S) in C&H,
                'Ks': 0.0077,   # saturated hydraulic conductivity, cm/min
                'S': 0.268      # sorptivity, cm/min^1/2  
            }
    
    Note: In C&H 1978, soil water retention relationships were defined according to _tensions_. 
    These tensions are specified as lengths, and are always _positive_ 
    (tension, like depth has an implied relationship to zero).

    To convert a tension, Psi_cm, (positive quantity of length) into a water potential, Psi_Pa, 
    (negative measure of energy density per unit volume, or Pa), you do the following:

    Psi_Pa = -1 * Psi_cm * rho * g

    This conversion is done during initiation of the soil class.

    """

    def __init__(self, texture=None, params=None):
        """ Initializes a soil object.

        The init function requires _either_ a soil texture or a params dictionary
        (see class description)
        
        """

        self._valid_params = set(['b', 'Psi_S_cm', 'Psi_l', 'n', 'Ks', 'S'])
        self._required_params = set(['b', 'Psi_S_cm', 'n', 'Ks'])
        
        # Set required attributes to None:
        [setattr(self, attr, None) for attr in self._required_params]
        
        if texture: # If this class is instanced with a specific USDA soil texture.
            texture = texture.lower()
            # Assign texture parameters based on the appropriate soil class:
            for attr, val in soils[texture].items():
                setattr(self, attr, val)

        elif params: # If the class is instanced with a set of soil parameters
            for attr, val in params.items():
                # Only include valid soil parameters
                if attr in self._valid_params:  
                    setattr(self, attr, val)

            # Check that all required parameters have been set
            if not self._required_params.issubset(self.__dict__.keys()):
                missing = self._required_params.difference(self.__dict__.keys())
                raise AttributeError("Missing required parameters, {list}".format(list=missing))
            
        else: 
            raise AttributeError("Must pass either a soil texture or dict of parameters")
        
        # Set Psi_S (MPa) from Psi_S_cm (cm). Assumes that Psi_S_cm is positive (as it should be!)
        self.Psi_S_MPa = -1 * self.Psi_S_cm / 100 * rho * g / 1E6 
        
        # Set Ks (mm/day) from Ks (cm/min).
        self.Ks = self.Ks*10*60*24

        # This version of sfc calculation comes from Laio et al. 2001b. Specifically, cf. the discussion
        # on p.714, and equation 15. 
        # self.sfc = pow(0.05/(self.Ks/10),1/(2*self.b+3))  # Convert Ks in cm/day 

        # This version of the sfc calculation uses the psi-theta relationships in Clapp & Hornberger to 
        # determine s_fc based on a texture-specific field_capacity.
        self.sfc = self.s(psi=field_capacity)
        
        # Make sure that field capacity is always lower than soil porosity.
        if self.sfc > 1:
            raise ValueError("soil field capacity, {sfc} is larger than 1".format(
                sfc=self.sfc
            ))

        # Set parameters related to pore size distribution index:
        self.Beta = 2*self.b + 4
        self.c = 2*self.b + 3

        # Hygroscopic point is when soil is so dry no further evaporation will occur.
        self.sh = self.s(self.theta(-10))               # Hygroscopic point in relative soil moisture [0-1]
        self.nZr = None                                 

    def _check_theta(self, theta):
        error = "theta, {theta}, must be be in the interval (0,{n}]".format(
                theta=theta, n=self.n)
        if theta > self.n or theta < 0:
            raise ValueError(error)

    def psi(self, theta):
        """ Return water potential in Pa based 
        on volumetric soil water content in m^3/m^3

        Note: Assumes that Psi is a water potential, and therefore Psi < 0 for unsaturated soils!

        Usage: psi(theta):
        
            theta = soil water content [m^3/m^3]
        
        """
        self._check_theta(theta)
        s = self.s(theta=theta)          
        return round(self.Psi_S_MPa * pow(s,-self.b),PRECISION)
    
    def theta(self,psi):
        """ Return a volumetric water content in m^3/m^3 
        based on a given water potential (MPa)

        Note: 
        Usage: theta(psi):

            psi = soil water potential [MPa]
        
        """
        if psi > 0:
            raise ValueError("psi, {psi}, must be less than or equal to zero.".format(psi=psi))
        # Ensure result is rounded to correct precision and that we do not exceed porosity
        return min([round((self.n * pow(psi/self.Psi_S_MPa, 1/-self.b)),PRECISION), self.n]) 

    def s(self,theta=None,psi=None):
        """ Return a relative soil moisture value, s [0-1]
        given a volumetric water content [m^3/m^3] or a 
        water potential [MPa]

        Usage: s(theta):

            theta = volumetric water content [m^3/m^3]
            psi = water potential [MPa]

        Note: theta must be in the interval 0-n (porosity)
        Note: psi must be negative
        Note: Function must be called with either theta or psi, but not both.
        
        """
        if theta and psi:
            raise ValueError(
            "Both theta ({theta}) and psi {psi} values provided only one argument allowed".format(
                theta=theta,
                psi=psi
            ))
        if psi:
            theta = self.theta(psi)
        self._check_theta(theta)
        try:
            return round(theta/self.n, PRECISION)
        except:
            raise ValueError("Either theta or psi must be provided as an argument.")

    def set_nZr(self,plant):
        # TODO: change this back to theta
        """ Sets the nZr for this soil in order to 
        determine fluxes in mm/day rather than relative
        soil moisture

        Usage: set_nZr(plant)

            plant = plant object with plant.Zr value set.
        
        Returns:

            nZr = n * Zr

            Also sets internal soil property nZr according to:

                self.nZr = self.n * plant.Zr
        """
        self.nZr = self.n * plant.Zr 
        return self.nZr 

    
    def calc_Q(self,s,units='mm/day'):
        """ Determines runoff as a function of relative soil moisture

        Usage: 

            calc_Q(s,units)

            s = relative soil moisture [0-1]
            units = units to return Drainage in
                options are 'mm/day' (default). 
                Otherwise, returns in [0-1] relative soil 
                moisture

        Returns:

            Q = runoff [mm/day] or [0-1]
        
        """

        # Initialize Q
        Q = 0
        
        # Saturation excess runoff occurs when 
        # relative soil moisture exceeds 1.
        if s > 1:
            Q = s - 1
        
        if units == 'mm/day':       
            return Q * self.nZr
        else:
            return Q

    def calc_D(self, s0, units='mm/day', timestep=1):
        """ Calculates daily drainage (Drainage) loss as a function of initial soil moisture in mm/day

        Usage: calc_D(s0, units)

            s0 = initial soil moisture [0-1]
            units = units to return Drainage in
                options are 'mm/day' (default). 
                Otherwise, returns in [0-1] relative soil 
                moisture

        Returns:

            D = Drainage [mm/day] if units='mm/day' 
            else returns Drainage in units of saturation [0-1]

        Notes: Drainage will be returned as a positive quantity.

        v1:
            - We are not using the calc_t_sfc code to avoid dependency on Emax
            - We use Dmax (s0-sfc) as an upper bound on the value of D that is returned.
            
        """
        
        Dmax = s0 - self.sfc # This is the largest amount of Drainage possible.  

        if Dmax > 0:
            # If the initial soil moisture, s0, is larger than field capacity, sfc,
            # Drainage is possible

            # Define m
            m =  self.Ks / (self.nZr * ( exp( self.Beta*(1 - self.sfc) ) - 1 ))  

            # Calculate D. 
            # Solution to D(s) = - ds/dt = - (s(t) - s(t0))
            # Note: this solution is correct. There is a missing close 
            # parenthesis in Laio et al. (2001) Eq. 20, so their equation should
            # not be used directly.
            D = (1 /  self.Beta ) * ln( 
                exp(self.Beta * (s0 - self.sfc)) 
                - exp(- m * self.Beta * timestep ) * (exp(self.Beta *  (s0 - self.sfc)) - 1) )

            D = min(D,Dmax) # Don't return a value larger than Dmax

            if units=='mm/day':
                return D * self.nZr
            else:
                # return in units of saturation/day
                return D
        else:
            # If the initial soil moisture, s0, is smaller than field capacity, sfc,
            # Drainage has ceased
            return 0



if __name__ == "__main__":
    import doctest
    doctest.testmod()  
# %%
