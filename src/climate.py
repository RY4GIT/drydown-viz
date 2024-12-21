# %%
"""
This code is based on the following source: 

natasha, Kelly Caylor, Noah Spahn, & Bryn Morgan. (2021). ecohydro/maize-Toff: First release of maize-Toff (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5204423

"""

#%% Climate Class Definition
import numpy as np
from numpy.random import exponential, uniform
from dateutil.relativedelta import *
from datetime import timedelta, datetime

class Climate():

    """ Creates a years worth of daily rainfall timeseries for use in ecohydrological modeling

    Usage: climate = Climate(alpha_r, lambda_r, ET_max)

        alpha_r = average storm depth [mm]
        lambda_r = storm frequency [day^-1]

    Default values:
        alpha_r = 10
        lambda_r = 0.25
        t_seas = 180
        ET_max = 6.5

    Note: lambda must either be a single value (constant rainfall probability all season),
    or have length of tseas (discrete rainfall probabilities each day.

    """
    def __init__(self, 
        alpha_r=[10.0] * 365, lambda_r=[0.25] * 365, lambda_std=[0.0]*365, 
        ET_max=6.5, q_e = 1.5, **kwargs):

        self.ET_max = ET_max
        self.q_e = q_e
        
        # Rainfall parameters
        self.alpha_r = alpha_r
        self.lambda_r = lambda_r
        self.lambda_std = lambda_std

        # Use the static method, generate, to create this instance's rainfall.
        self.rainfall = self.generate(self.alpha_r, self.lambda_r)


    def calc_E(self, s, LAI=None, sh=None): 
        """ Determines the daily evaporation as a function of relative soil moisture

        Usage: calc_E(s)

            s = relative soil moisture [0-1]

        """
        if LAI == None:
            raise ValueError("Climate calc_E expects LAI that's not None.")

        k = 0.5
        E_max_p = self.ET_max*np.exp(-k*LAI) 
        if s >= sh:
            return pow((s-sh)/(1-sh), self.q_e)*E_max_p
        else:
            return 0
            
    @staticmethod # Static methods can be called without instancing the class.
    def generate(alpha_r, lambda_r, t_sim=365, doy_start=1):
        """ Makes a time series of rainfall based on parameters

        Usage:
            generate(alpha_r, lambda_r, t_seas)

            alpha_r = average storm depth [mm]
            lambda_r = storm frequency [day^-1]
            t_sim = length of growing season [days]
            doy_start = day of year to start the simulation [day] 

        Note: lambda must either be a single value (constant rainfall probability all season),
        or have length of tseas (discrete rainfall probabilities each day. 

        """
        # Force doy to be in [1,365]:
        doys = np.arange(doy_start, doy_start + t_sim)
        while (doys - 365 > 0).any() == True:
            doys = doys - 365 * ((doys - 365) > 0)
        amounts = [exponential(scale=alpha_r[doy-1], size=1)[0] for doy in doys]
        rain_days = [(uniform(low=0, high=1, size=1) <= lambda_r[doy-1] ).astype(int) for doy in doys]
        return np.multiply(amounts, [v[0] for v in rain_days])