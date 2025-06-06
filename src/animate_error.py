# %% Animate Loss function vs Drydown curves in drainage stage (D)

__author__ = "Ryoko Araki"
__contact__ = "raraki8159@sdsu.edu"
__copyright__ = "Ryoko Araki"

__license__ = "MIT"
__status__ = "development"

# IMPORTS
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from functions import (
    theta_3stage,
    lossfunc_3stage,
    loss_drainage,
    loss_ET_stageI,
    theta_drainage,
    find_t_fc,
)

# Math font
import matplotlib as mpl

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
mpl.rcParams["font.size"] = 14.0
mpl.rcParams["axes.titlesize"] = 12.0
plt.rcParams["mathtext.fontset"] = "stixsans"

# Define plot labels
var_dict = {
    "theta": {
        "symbol": r"$\theta$",
        "label": r"Soil moisture, $\theta$",
        "unit": r"($\mathrm{m}^3$ $\mathrm{m}^{-3}$)",
    },
    "dtheta": {
        "symbol": r"$-\frac{d\theta}{dt}$",
        "label": r"$\minus\frac{d\theta}{dt}$",
        "unit": r"($\mathrm{m}^3$ $\mathrm{m}^{-3}$ $\mathrm{day}^{-1}$)",
    },
    "t": {
        "symbol": r"$t$",
        "label": r"$t$",
        "unit": r"(day)",
    },
}
theta_vardict = var_dict["theta"]
dtheta_vardict = var_dict["dtheta"]
# %% Define parameters

# ______________________
# Soil hydraulic parameters
Ks = 30  # Saturated hydraulic conductivity in mm/day
n = 0.6  # Porosity in m3/m3.
beta = 32  # Parameter b
z = 50  # Soil thickness in mm (default is 50 mm).


theta_fc = 0.40  # Field capacity in m3/m3
theta_star = 0.35  # Critical soil moisture content in m3/m3
theta_w = 0.05  # Wilting point soil moisture content in m3/m3

# ______________________
# Climate & plant parmaeters
q = 1  # Degree of non-linearity in the soil moisture response related to plant stress
ETmax = 6  # Maximum evapotranpisration rate in mm/day.
k = ETmax / z  # Normalized ETmax rate per soil depth

# %%
# Define variables

# _________________________________________________________________________________
# VARIABLES FOR BASE PLOTS (entire shape of the loss function and drydowns)
delta = 1e-03
theta_base = np.arange(0, n, delta)
theta_0_base = 0.6
tmax_base = 10
t_base = np.arange(0, tmax_base, delta)
# Drainage stage
# theta_D = np.arange(theta_fc, n, 1e-03)

loss_values_base = [
    lossfunc_3stage(
        theta=theta_t,
        q=q,
        ETmax=ETmax,
        Ks=Ks,
        beta=beta,
        n=n,
        theta_fc=theta_fc,
        theta_star=theta_star,
        theta_w=theta_w,
        z=z,
    )
    for theta_t in theta_base
]

t_fc_base = find_t_fc(
    ETmax=ETmax, Ks=Ks, beta=beta, theta_fc=theta_fc, theta_0=theta_0_base, n=n, z=z
)

dd_values_base = [
    theta_3stage(
        t=t_i,
        q=q,
        ETmax=ETmax,
        theta_0=theta_0_base,
        theta_w=theta_w,
        theta_star=theta_star,
        theta_fc=theta_fc,
        Ks=Ks,
        beta=beta,
        n=n,
    )
    for t_i in t_base
]

# %%
# _________________________________________________________________________________
# VARIABLES FOR SAMPLE PLOTS

t_base_sample = np.arange(0, tmax_base, 1.0)

_dd_values_sample = [
    theta_3stage(
        t=t_i,
        q=q,
        ETmax=ETmax,
        theta_0=theta_0_base,
        theta_w=theta_w,
        theta_star=theta_star,
        theta_fc=theta_fc,
        Ks=Ks,
        beta=beta,
        n=n,
    )
    for t_i in t_base_sample
]

# Add Gaussian noise

# Parameters for Gaussian noise
mean = 0  # Mean of the noise
std_dev = 0.03  # Standard deviation of the noise
np.random.seed(42)

noise = np.random.normal(mean, std_dev, len(_dd_values_sample))
dd_values_sample = _dd_values_sample + noise

loss_values_sample = np.diff(dd_values_sample)
# %%
# _________________________________________________________________________________
# Set up figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
dtheta_color = "#FF8383"
base_color = "#D3D3D3"
linewidth = 3


# Animation update function
def update(frame):
    #################################################################################
    # Preparation
    #################################################################################

    # _________________________________________________________________________________
    # Check exit conditions
    if frame > len(dd_values_base) - 1:
        return []

    #################################################################################
    # Update Loss function plot
    #################################################################################
    ax1.clear()

    # _________________________________________________________________________________
    # base plot
    ax1.plot(
        theta_base,
        loss_values_base,
        linewidth=linewidth,
        color=base_color,
        alpha=0.5,
        label="Theoretical",
    )

    # _________________________________________________________________________________
    # Error plot
    ax1.scatter(
        dd_values_sample[:-1][:frame],
        loss_values_sample[:frame],
        marker="o",
        color=dtheta_color,
        label="Sample + noise",
    )

    # _________________________________________________________________________________
    # Adjustment
    ax1.invert_yaxis()
    ax1.set_xlim([0.0, n])
    ax1.set_ylim([0.0, -Ks / z])
    ax1.set(
        xlabel=f"{theta_vardict['label']} {theta_vardict['unit']}",
        ylabel=f"{var_dict['dtheta']['label']} {var_dict['dtheta']['unit']}",
    )
    ax1.set_title("Loss function space", loc="left")
    ax1.set_xticks(
        [theta_w, theta_star, theta_fc, n],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$", r"$\theta_{\mathrm{fc}}$", r"$n$"],
    )
    ax1.set_yticks(
        [0, -k, -Ks / z],
        [0, r"$\frac{\mathrm{ET}_{\mathrm{max}}}{z}$", r"$\frac{Ks}{z}$"],
    )
    ax1.legend(frameon=False)

    #################################################################################
    # Update Drydown plot
    #################################################################################

    ax2.clear()

    # _________________________________________________________________________________
    # base plot
    ax2.plot(
        t_base,
        dd_values_base,
        linewidth=linewidth,
        color=base_color,
        alpha=0.5,
        label="Theoretical",
    )

    # _________________________________________________________________________________
    # Error plot
    ax2.scatter(
        t_base_sample,
        dd_values_sample,
        marker="o",
        color=base_color,
    )

    ax2.plot(
        t_base_sample[:frame],
        dd_values_sample[:frame],
        linewidth=linewidth / 2,
        color=dtheta_color,
        alpha=0.5,
        label="Sample + noise",
    )

    # _________________________________________________________________________________
    # Adjustment
    ax2.set_ylim([0.0, n])
    ax2.set_xlim([0, tmax_base])
    ax2.set(
        xlabel=f"{var_dict['t']['label']} {var_dict['t']['unit']}",
        ylabel=f"{theta_vardict['label']} {theta_vardict['unit']}",
    )
    ax2.set_yticks(
        [theta_w, theta_star, theta_fc, n],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$", r"$\theta_{\mathrm{fc}}$", r"$n$"],
    )
    ax2.set_title("Drydown space", loc="left")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(frameon=False)

    fig.tight_layout()

    return ax1, ax2


# _________________________________________________________________________________
# Create animation
ani = FuncAnimation(fig, func=update, frames=12, interval=300, blit=False)
# frames: The total number of frames for the animation.
# interval: The time between frames in milliseconds(e.g., interval=50 for 20 FPS)

# _________________________________________________________________________________
# Save animation
out_dir = "./out"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
out_filename = "soil_est_ver1.gif"
ani.save(os.path.join(out_dir, out_filename), writer="pillow", dpi=300)
print(f"Drainage stage animation was saved to: {os.path.join(out_dir, out_filename)}")

# %%
