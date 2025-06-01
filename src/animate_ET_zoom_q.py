# %% Animate Loss function vs Drydown curves in drainage stage (D)

__author__ = "Ryoko Araki"
__contact__ = "raraki8159@sdsu.edu"
__copyright__ = "Ryoko Araki"

__license__ = "MIT"
__status__ = "development"

# %% IMPORTS
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from functions import (
    theta_3stage,
    lossfunc_3stage,
    loss_ET_stageI,
    loss_ET_stageII,
    theta_ET_stageI,
    theta_ET_stageII,
    find_t_fc,
    find_t_star,
)

# Math font
import matplotlib as mpl

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
mpl.rcParams["font.size"] = 12.0
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


q_agg = 0.7
q_csv = 1.5
# %%
# Define variables

# _________________________________________________________________________________
# VARIABLES FOR BASE PLOTS (entire shape of the loss function and drydowns)
delta = 1e-03
theta_base = np.arange(0, theta_fc, delta)
theta_0_base = 0.6
tmax_base = 15
t_base = np.arange(0, tmax_base, delta)

# Drainage stage
# theta_ET = np.arange(theta_fc, n, 1e-03)

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
t_star_base = (
    find_t_star(theta_0=theta_fc, ETmax=ETmax, theta_star=theta_star) + t_fc_base
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
# Set up figure
fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(9, 4))
ET_color = "#5DB996"
ET_color_agg = "#F7CA0D"
ET_color_csv = "#01665E"
base_color = "#D3D3D3"
linewidth = 3


# Animation update function
def update(frame):
    #################################################################################
    # Preparation
    #################################################################################

    # _________________________________________________________________________________
    # Define the timerange
    global t_Delta
    t_Delta = (frame + 1) * 0.02  # Increment for each frame
    t_ET = np.arange(t_star_base, t_star_base + t_Delta, 1e-2)

    # _________________________________________________________________________________
    # DRAINAGE STAGE PLOT
    dd_values_ET = theta_ET_stageII(
        t=t_ET,
        q=q,
        ETmax=ETmax,
        theta_0=theta_star,
        theta_star=theta_star,
        theta_w=theta_w,
        t_star=t_star_base,
        z=z,
    )

    dd_values_ET_agg = theta_ET_stageII(
        t=t_ET,
        q=q_agg,
        ETmax=ETmax,
        theta_0=theta_star,
        theta_star=theta_star,
        theta_w=theta_w,
        t_star=t_star_base,
        z=z,
    )

    dd_values_ET_csv = theta_ET_stageII(
        t=t_ET,
        q=q_csv,
        ETmax=ETmax,
        theta_0=theta_star,
        theta_star=theta_star,
        theta_w=theta_w,
        t_star=t_star_base,
        z=z,
    )

    # Check exit conditions
    if dd_values_ET[-1] < theta_w:
        return []

    theta_ET = np.arange(dd_values_ET[-1], theta_star, 1e-03)

    loss_values_ET = loss_ET_stageII(
        theta=theta_ET,
        ETmax=ETmax,
        theta_star=theta_star,
        theta_w=theta_w,
        q=1.0,
        z=50.0,
    )

    loss_values_ET_agg = loss_ET_stageII(
        theta=theta_ET,
        ETmax=ETmax,
        theta_star=theta_star,
        theta_w=theta_w,
        q=q_agg,
    )

    loss_values_ET_csv = loss_ET_stageII(
        theta=theta_ET,
        ETmax=ETmax,
        theta_star=theta_star,
        theta_w=theta_w,
        q=q_csv,
    )

    #################################################################################
    # Update Loss function plot
    #################################################################################
    ax1.clear()

    # _________________________________________________________________________________
    # base plot
    ax1.plot(
        theta_base, loss_values_base, linewidth=linewidth, color=base_color, alpha=0.5
    )

    # _________________________________________________________________________________
    # Drainage plot
    ax1.plot(theta_ET, loss_values_ET, linewidth=linewidth, color=ET_color)
    ax1.plot(
        theta_ET, loss_values_ET_agg, linewidth=linewidth, color=ET_color_agg, alpha=0.5
    )
    ax1.plot(
        theta_ET, loss_values_ET_csv, linewidth=linewidth, color=ET_color_csv, alpha=0.5
    )

    # ax1.scatter(theta_ET[0], loss_values_ET[0], color=ET_color)
    # ax1.scatter(theta_ET[0], loss_values_ET_agg[0], color=ET_color, alpha=0.5)
    # ax1.scatter(theta_ET[0], loss_values_ET_csv[0], color=ET_color, alpha=0.5)

    # _________________________________________________________________________________
    # Adjustment
    ax1.invert_yaxis()
    ax1.set_xlim([0.0, theta_fc * 1.1])
    ax1.set_ylim([0.0, -k * 1.1])
    ax1.set(
        xlabel=f"{theta_vardict['label']} {theta_vardict['unit']}",
        ylabel=f"{var_dict['dtheta']['label']} {var_dict['dtheta']['unit']}",
    )
    # ax1.set_title("Loss function", loc="left", fontweight="bold")
    ax1.set_xticks(
        [theta_w, theta_star, theta_fc],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$", r"$\theta_{\mathrm{fc}}$"],
    )
    ax1.set_yticks(
        [0, -k],
        [0, r"$\frac{\mathrm{ET}_{\mathrm{max}}}{z}$"],
    )

    #################################################################################
    # Update Drydown plot
    #################################################################################

    ax2.clear()

    # _________________________________________________________________________________
    # base plot
    ax2.plot(t_base, dd_values_base, linewidth=linewidth, color=base_color, alpha=0.5)

    # _________________________________________________________________________________
    # Drainage plot
    ax2.plot(t_ET, dd_values_ET, linewidth=linewidth, color=ET_color)
    # ax2.scatter(t_ET[-1], dd_values_ET[-1], color=ET_color)
    ax2.plot(t_ET, dd_values_ET_agg, linewidth=linewidth, color=ET_color_agg, alpha=0.5)
    # ax2.scatter(t_ET[-1], dd_values_ET_agg[-1], color=ET_color, alpha=0.5)
    ax2.plot(t_ET, dd_values_ET_csv, linewidth=linewidth, color=ET_color_csv, alpha=0.5)
    # ax2.scatter(t_ET[-1], dd_values_ET_csv[-1], color=ET_color, alpha=0.5)

    # _________________________________________________________________________________
    # Adjustment
    ax2.set_ylim([0.0, theta_fc * 1.1])
    ax2.set_xlim([0, tmax_base])
    ax2.set_xticks(np.arange(0, tmax_base + 1, 5))
    ax2.set(
        xlabel=f"{var_dict['t']['label']} {var_dict['t']['unit']}",
        ylabel=f"{theta_vardict['label']} {theta_vardict['unit']}",
    )
    ax2.set_yticks(
        [theta_w, theta_star, theta_fc],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$", r"$\theta_{\mathrm{fc}}$"],
    )
    # ax2.set_title("Drydown", loc="left", fontweight="bold")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()

    return ax1, ax2


# _________________________________________________________________________________
# Create animation
ani = FuncAnimation(fig, update, frames=600, interval=25, blit=False)
# frames: The total number of frames for the animation.
# interval: The time between frames in milliseconds(e.g., interval=50 for 20 FPS)

# _________________________________________________________________________________
# Save animation
out_ETir = "./out"
if not os.path.exists(out_ETir):
    os.mkdir(out_ETir)
out_filename = "soil_ET_zoom_q.gif"
ani.save(os.path.join(out_ETir, out_filename), writer="pillow", dpi=300)
print(f"Drainage stage animation was saved to: {os.path.join(out_ETir, out_filename)}")

# %%
