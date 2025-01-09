# %%
"""
This code is based on the following source: 

natasha, Kelly Caylor, Noah Spahn, & Bryn Morgan. (2021). ecohydro/maize-Toff: First release of maize-Toff (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5204423

"""

# %% Plant Class Definition
import numpy as np


class Plant():
    """ Defines a plant class


    """

    def __init__(self,
                 Zr=400,             # Rooting depth [mm]
                 T_max=4,            # Maximum transpiration [mm/day]
                 # wilting point of plant in water potential [MPa]
                 sw_MPa=-1.5,
                 # water potential of maximum transpiration [MPa],
                 s_star_MPa=-0.05,
                 soil=None          # a soil in which this plant will grow
                 ):
        self.Zr = Zr
        self.sw_MPa = sw_MPa
        self.s_star_MPa = s_star_MPa

        # Use the soil to determine critical plant parameters
        self.sw = soil.s(soil.theta(sw_MPa))
        self.s_star = soil.s(soil.theta(s_star_MPa))

        self.T_max = T_max      # Should be over-written by any subclass.

    def calc_LAI(self):
        raise NotImplementedError

    def calc_T(self, s):
        raise NotImplementedError

    def calc_I(self):
        raise NotImplementedError

# %%


class Crop(Plant):
    """ Creates a Crop class.

    Usage: crop = Crop(soil=soil)

    optional keyword arguments and their default values:
        'Zr': 500,          # Planting depth [mm]
        'sw_MPa':-1.5,      # Plant wilting point [MPa]
        's_star_MPa':-0.2,  # Water potential for max T
        'kc_max':1.2,       # Maximum crop coefficient
        'LAI_max':3.0,      # Max Leaf Area Index [m2/m2]
        'T_max':4.0         # Max Crop Water Use [mm/day]

    """

    def __init__(self, kc_max=1.2, LAI_max=3.0, T_max=4, lgp=180, F1=0.16, F2=0.44, F3=0.76,
                 eos=1.0, kc_ini=0.30, kc_eos=0.6, q=1, *args, **kwargs):

        # Maximum crop coefficient; Kc at Reproductive Stage [0-1]
        self.kc_max = kc_max
        self.LAI_max = LAI_max    # Maximum crop leaf area index [m^2/m^2]
        self.T_max = T_max        # Maximum crop water use [mm/day]
        self.lgp = lgp            # Length of growing period [days]
        self.F1 = F1              # Fraction of Season from Initial to Vegetative
        self.F2 = F2              # Fraction of Season from Initial to Reproductive
        self.F3 = F3              # Fraction of Season from Initial to Ripening
        self.eos = eos            # Fraction of Season at End
        self.kc_ini = kc_ini      # Kc at Initial Stage
        self.kc_eos = kc_eos      # Kc at End of Season
        self.q = q                # Q_stress for static stress calculation
        super(Crop, self).__init__(*args, **kwargs)

    def calc_kc(self, day_of_season=0):
        """ Calculates crop coefficient that varies throughout the season 

        Usage: calc_kc(self, day_of_season)

            day_of_season = 0 # Start date [day]

        Note: t must be a single-dimension array. day_of_season is when the plant starts
        to grow and is hard coded to zero because python is zero-indexed.
        Kc values come from Table 11. Lengths of crop development stages for 
        maize (grain) in East Africa (altitude), Chapter 6 in FAO (1998).

        """
        if not day_of_season >= 0:
            raise ValueError("day_of_season must be >= 0")
        if day_of_season <= self.lgp*self.F1:
            return self.kc_ini
        elif day_of_season < self.lgp*self.F2:
            return ((self.kc_max-self.kc_ini)/(self.F2*self.lgp-self.F1*self.lgp))*(day_of_season-self.F1*self.lgp)+self.kc_ini
        elif day_of_season <= self.lgp*self.F3:
            return self.kc_max
        elif day_of_season < self.lgp*self.eos:
            return ((self.kc_eos-self.kc_max)/(self.eos*self.lgp-self.F3*self.lgp))*(day_of_season-self.F3*self.lgp)+self.kc_max
        else:
            return self.kc_eos

    def calc_T_max(self, kc):
        """ Calculates max Transpiration variable.

        Usage: calc_T_max(kc)

        T = kc * T_max

        """
        return kc * self.T_max

    def _kc_from_LAI(self, LAI, p=1):
        """ Returns a kc variable. kc comes
            from function of LAI. Currently based on linear relationship 
            between kc and LAI (assumption).

        Usage: _kc_from_LAI(LAI, p=1)

            kc = (kc_max/LAI_max)^p * LAI.

        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return np.power((self.kc_max/self.LAI_max), p) * LAI

    def calc_LAI(self, kc, p=1):
        """ Returns a Leaf Area Index (LAI) variable. LAI comes
            from function of kc. Currently based on linear relationship 
            between kc and LAI (assumption).

        Usage: calc_LAI(kc, p=1)

            LAI = (LAI_max/kc_max)^p * kc.

        Note: p=1 assumes a linear relationship between LAI and kc

        """
        return np.power((self.LAI_max/self.kc_max), p) * kc

    def calc_T(self, s, LAI=None, kc=None):
        """ Calculates Transpiration variable as a stepwise
            linear function. Will use LAI value if both LAI 
            and kc are provided.

        Usage: calc_T(s, LAI, kc)

            s = relative soil moisture [0-1]
            LAI = leaf area index [m2/m2]
            kc = crop coefficient [-].

        Note: Either LAI or kc must be provided.

        """
        # if not LAI and not kc:
        # raise(ValueError, "Function requires either LAI or kc to be set.")
        if LAI:
            kc = self._kc_from_LAI(LAI)
        if kc:
            if s >= self.s_star:
                return self.calc_T_max(kc)
            elif s >= self.sw:
                return self.calc_T_max(kc) * np.power(((s-self.sw)/(self.s_star-self.sw)), self.q)
            else:
                return 0

    def calc_I(self, LAI, int_efficiency=1):
        """
        Determines canopy interception based on crop LAI and
        an interception efficiency term.

        Usage: calc_I(LAI, int_efficiency)

            LAI = leaf area index [m2/m2]
            int_efficiency = conversion term [-].

        """
        return LAI * int_efficiency  # 0 for no int

    def calc_stress(self, s):
        """ Calculates static water stress.

        Usage: calc_stress(s, q=2)

            s = relative saturation [0-1]
            q = 2

        Note: The value of q changes based on plant species or soil type. 
        See equation 4.13, p.101 in Rodriguez-Iturbe & Porporato (2004)
        """
        if s < self.sw:
            stress = 1
        elif s >= self.s_star:
            stress = 0
        else:
            stress = ((self.s_star - s)/(self.s_star - self.sw))**self.q
        return stress

    # Y_MAX was previously set to 4260 for 180 day crop
    def calc_dstress(self, s, stress, Y_MAX=None):
        '''Calculates dyamic water stress (theta) which is a measure of total water stress during the growing season
        as proposed in Porporato et al. (2001). Considers the duration and frequency of water defict periods below a 
        critical value. The function also calculates yield based on dynamic water stress and returns three items in a
        list: average static water stress, dynamic water stress, and yield in kg per ha. 

        Usage: calc_dstress(s, stress):

            s = relative saturation [0-1]
            stress = static stress [0-1]
            Y_MAX = None [t / ha]

        Default values:
            mstr_memb = average static stress [0-1]
            mcrs_memb = average duration of water stress [days]
            ncrs_memb = average frequency of water stress [dim]
            self.lgp = length of the growing season [days]
            K_PAR = fraction of growing season before crop fails [dim]
            R_PAR = effect of number of excursions below stress point [dim]
            INVL_SIMU = number of daily timesteps used in calculating the soil moisture time series [dim]

        Returns:
            mstr_memb = np.mean(((self.s_star - s)/(self.s_star - self.sw))**q) # average static water stress
            dstr_memb = (mstr_memb * mcrs_memb) / (K_PAR * self.lgp))**(ncrs_memb**-R_PAR) # dynamic water stress
            yield_kg_ha = Y_MAX * (1 - dstr_memb) # yield in kg per ha

        Note: Y_max can be set to the potential yield of the variety of interest. During implementation we use the 
        Y_max set by the linear regression based on the data from Kenya Seed Co. with this line of code
        data = [crop.calc_dstress(s=df.s, stress=df.stress, Y_MAX = evolved_calc_yield(dtm=lgp)) for df in output]

        '''
        # if Y_MAX = NotImplementedError:
        #        raise ValueError("lambda_r values and alpha_r values must be same length")

        # Step 0. Define variables
        # K parameter is the portion of the season that crop can endure stress before it fails.
        K_PAR = 0.25
        # R parameter should not change since it is the square root term in Porporato et al. (2001) p. 739
        R_PAR = 0.5
        INVL_SIMU = 1

        # Step 1. Calculate average static stress
        if len(stress) > 0:
            # Subset the growing period and get avg soil moisture
            start = 60
            end = start + self.lgp
            stress_subset = stress[start:end]
            mstr_memb = np.mean(stress_subset)
        else:
            mstr_memb = 0.

        # Step 2. Calculate threshold crossing parameters
        # Select indices of s time series where s is below wilting point
        indx_memb = np.where(s >= self.s_star)
        # Append to an array using np.append where last value is lgp+1 and INVL_SIMU is how many simulations are being run
        # Then have zero be the first item and
        # with np.diff give the difference to find the soil moisture difference between s_star and the excursion
        # play around with this to figure it out
        ccrs_memb = np.diff(np.append(0, np.append(
            indx_memb, INVL_SIMU * self.lgp + 1))) - 1
        # The duration of water stress events where there is stress because value is greater than 0
        ccrs_memb = ccrs_memb[ccrs_memb > 0]
        # Variable with number of excursions below wilting point (frequency)
        ncrs_memb = len(ccrs_memb)  # dim
        if ncrs_memb > 0:
            # if there are more than 0 excursions then calculate mean of duration of water stress and divide by INVL_SIMU
            mcrs_memb = np.mean(ccrs_memb) / INVL_SIMU  # days
        else:
            mcrs_memb = 0.

        # Step 3. Calculate dynamic stress
        dstr_memb = ((mstr_memb * mcrs_memb) /
                     (K_PAR * self.lgp))**(ncrs_memb**-R_PAR)
        if dstr_memb > 1.:
            dstr_memb = 1.

        # Step 4. Calculate yield
        yield_kg_ha = Y_MAX * (1 - dstr_memb)

        return mstr_memb, dstr_memb, yield_kg_ha


class StaticCrop(Crop):
    """ Creates a StaticCrop class that inherits from Crop and overrides the calc_kc method.

    Usage: crop = StaticCrop(soil=soil)

    """

    def __init__(self, const_kc=0.7, q_t=1, *args, **kwargs):
        self.const_kc = const_kc
        self.q_t = q_t
        super(StaticCrop, self).__init__(*args, **kwargs)

    def calc_kc(self, day_of_season):
        """ 
        Calculates crop coefficient that does NOT vary throughout the season
        """
        return self.const_kc

    def calc_T(self, s, LAI=None, kc=None):
        """ Calculates Transpiration variable as a stepwise
            linear function. Will use LAI value if both LAI 
            and kc are provided.

        Usage: calc_T(s, LAI, kc)

            s = relative soil moisture [0-1]
            LAI = leaf area index [m2/m2]
            kc = crop coefficient [-].

        Note: Either LAI or kc must be provided.

        """
        # if not LAI and not kc:
        # raise(ValueError, "Function requires either LAI or kc to be set.")
        if LAI:
            kc = self._kc_from_LAI(LAI)
        if kc:
            if s >= self.s_star:
                return self.calc_T_max(kc)
            elif s >= self.sw:
                return np.power((s-self.sw)/(self.s_star-self.sw), self.q_t)*self.calc_T_max(kc)
            else:
                return 0
