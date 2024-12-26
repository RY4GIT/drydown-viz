# https://github.com/RY4GIT/maize-Toff/blob/master/notebooks/tau_playground3_simplerT.ipynb

# %%
import numpy as np
from model import ZeroECropModel, CropModel
from climate import Climate
from soil import Soil
from plant import Plant, Crop, StaticCrop
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rcParams

# %% ====================================
# MODEL SETTING
# =======================================

texture = "clay loam"
alpha_r = 10.0
lambda_r = 0.25
lambda_std = 0.0
ETmax = 6.5
q_e = 1.5

s0=0.3
planting_date=100
t_before=60
t_after=60
t_warmup=50

# %% ====================================
# MODEL RUN
# =======================================
soil = Soil(texture=texture)
climate= Climate(
        alpha_r=[alpha_r] * 365, 
        lambda_r=[lambda_r] * 365, 
        lambda_std=[lambda_std]*365, 
        ET_max=ETmax, 
        q_e = q_e
        )
plant = StaticCrop(soil=soil, q_t=1.5)

cropmodel = ZeroECropModel(soil=soil,
            climate=climate,
            crop=plant
            )

cropmodel.run(
            s0=s0,
            planting_date=planting_date,
            t_before=t_before,
            t_after=t_after,
            t_warmup=t_warmup
            )

# %% ====================================
# VISUALIZE
# =======================================

## Figure to animate: Time series of rainfall, soil moisture, and loss function 
plt.rcParams["font.family"] = "DejaVu Sans"  # Or any other available font
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]  # Ensure the font is set correctly
# mpl.rcParams["font.family"] = "sans-serif"
# mpl.rcParams["font.sans-serif"] = "Myriad Pro"
rcParams["font.size"] = 16.0
rcParams["axes.titlesize"] = 16.0
plt.rcParams["mathtext.fontset"] = (
    "stixsans"  #'stix'  # Or 'cm' (Computer Modern), 'stixsans', etc.
)

# =======================================
# Create a figure with GridSpec layout
fig = plt.figure(figsize=(18, 6)) 
gs = GridSpec(2, 6, figure=fig, width_ratios=[1, 1, 1, 1, 1, 4])  

# Create subplots for the large plot (2 rows, 1 column)
ax1 = fig.add_subplot(gs[0, 0:5])  # Top subplot for rainfall
ax2 = fig.add_subplot(gs[1, 0:5])  # Bottom subplot for soil moisture

# Create the small square plot on the right
ax3 = fig.add_subplot(gs[:, 5])  # Small plot spans both rows in the last column

# =======================================
# Figure a: Rainfall timeseries from simulation
output = cropmodel.output()
rf = output['R']
stress = output['stress']
rf.index = rf.index+1 # zero-indexed --> day of the season 
stress.index = stress.index+1

ax1.bar(rf.index, rf, color='k', edgecolor='k') 
ax1.set_ylabel='Daily rainfall [mm]'
ax1.set_xlabel='Day of season [day]'

# =======================================
# Figure b: Soil moisture timeseries from simulation
s = output['s']
s_sorted = np.sort(s)
_sstar = plant.s_star
_sw = plant.sw
s.index = s.index + 1 # zero-indexed --> day of the season 

porosity = Soil(texture).n 
theta = s * porosity
theta_sorted= np.sort(theta)
_theta_star = plant.s_star*porosity
_theta_sw = plant.sw*porosity

# LEGEND
hs = []
cycle = ['-','--']
var = [_sstar, _sw]
LABEL = ['Stress point (s*)','Wilting point (sw)']


ax2.plot(theta, color='k', lw=1)
ax2.set_ylabel(r'Soil moisture, $\theta$')
ax2.set_xlabel('Day of season [day]')

# =======================================
# Figure c: Soil moisture loss function

# ax3.plot(theta_sorted, [plant.calc_T(x, LAI=1.5) for x in s_sorted], '-', linewidth=2)
ax3.plot(theta_sorted, [soil.calc_Q(x)+soil.calc_D(x)+plant.calc_T(x, LAI=1.5) for x in s_sorted], '-', linewidth=2)
ax3.set_xlabel(r'Soil moisture, $\theta$') 
ax3.set_ylabel(r'Loss rate, $-d\theta/dt$ [mm/day]') 
ax3.set_ylim(-0.1,5.5)
ax3.set_xlim(_theta_sw-0.05,porosity+0.05)

fig.tight_layout()


# %%
output
# %%
