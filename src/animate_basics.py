# %%
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
    find_t_fc
)

# Math font
import matplotlib as mpl
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
mpl.rcParams["font.size"] = 14.0
mpl.rcParams["axes.titlesize"] = 12.0
plt.rcParams["mathtext.fontset"] = ("stixsans")
# %%
# Define parameters
q = 1

ETmax = 6
z = 50
k = ETmax / z

Ks = 30  # mm/day # Hydraulic conductivity constant
n = 0.6
beta = 32
theta_fc = 0.40
theta_star = 0.35
theta_w = 0.05


# Define plot format
var_dict = {
    "theta": {
        "symbol": r"$\theta$",
        "label": r"Soil moisture, $\theta$",
        "unit": r"($\mathrm{m}^3$ $\mathrm{m}^{-3}$)",
        "lim": [0, n],
    },
    "dtheta": {
        "symbol": r"$-\frac{d\theta}{dt}$",
        "label": r"$\minus\frac{d\theta}{dt}$",
        "unit": r"($\mathrm{m}^3$ $\mathrm{m}^{-3}$ $\mathrm{day}^{-1}$)",
        "lim": [-0.08, 0],
    },
    "t": {
        "symbol": r"$t$",
        "label": r"$t$",
        "unit": r"(day)",
    },
}
theta_vardict = var_dict["theta"]
dtheta_vardict = var_dict["dtheta"]

plt.rcParams["mathtext.fontset"] = (
    "stixsans"  # 'stix'  # Or 'cm' (Computer Modern), 'stixsans', etc.
)

# %%

# %%
# Define variables
theta = np.arange(0, n, 1e-03)
theta_0 = 0.6

tmax = 15
t_star = (theta_fc - theta_star) / k
t_before_star = np.arange(0, t_star, 1e-03)
t = np.arange(0, tmax, 1e-03)

# Drainage stage
# theta_D = np.arange(theta_fc, n, 1e-03)

loss_values = [lossfunc_3stage(
    theta=theta_t, q=q, ETmax=ETmax, Ks=Ks, beta=beta, n=n, theta_fc=theta_fc, theta_star=theta_star, theta_w=theta_w, z=z) for theta_t in theta]

# loss_values_D = [loss_drainage(
# theta=theta_t, Ks=Ks, beta=beta, theta_fc=theta_fc, n=n, z=50) + loss_ET_stageI(ETmax=ETmax, z=z) for theta_t in theta_D]

t_fc = find_t_fc(
    ETmax=ETmax, Ks=Ks, beta=beta, theta_fc=theta_fc, theta_0=theta_0, n=n, z=z)

# t_D = np.arange(0, t_fc, 1e-02)

dd_values = [theta_3stage(
    t=t_i, q=q, ETmax=ETmax, theta_0=theta_0, theta_w=theta_w, theta_star=theta_star, theta_fc=theta_fc, Ks=Ks, beta=beta, n=n) for t_i in t]

# dd_values_D = [theta_drainage(
# t=t_i, ETmax=ETmax, Ks=Ks, beta=beta, theta_fc=theta_fc, theta_0=theta_0, n=n) for t_i in t_D]


# %%
# # Set up figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
D_color = "#81BFDA"
base_color = "#D3D3D3"
linewidth = 3


# Animation update function


def update(frame):
    global Delta
    Delta = (frame + 1) * 0.02  # Increment for each frame
    t_D = np.arange(0, Delta, 1e-2)

    # Check exit conditions
    if (Delta > t_fc):
        # ani.event_source.stop()  # Stop the animation gracefully
        return []

    dd_values_D = [
        theta_drainage(t=t_i, ETmax=ETmax, Ks=Ks, beta=beta,
                       theta_fc=theta_fc, theta_0=theta_0, n=n)
        for t_i in t_D
    ]

    theta_D = np.arange(dd_values_D[-1], n, 1e-03)

    loss_values_D = [
        loss_drainage(theta=theta_t, Ks=Ks, beta=beta,
                      theta_fc=theta_fc, n=n, z=z)
        + loss_ET_stageI(ETmax=ETmax, z=z)
        for theta_t in theta_D
    ]

    # Update Loss function plot
    ax1.clear()
    # base plot
    ax1.plot(
        theta,
        loss_values,
        linewidth=linewidth,
        color=base_color,
        alpha=0.5
    )
    # Drainage plot
    ax1.plot(theta_D, loss_values_D, linewidth=linewidth, color=D_color)
    ax1.scatter(theta_D[0], loss_values_D[0], color=D_color)
    ax1.invert_yaxis()
    ax1.set_xlim([0.0, n])
    ax1.set_ylim([0.0, -Ks/z])
    ax1.set(
        xlabel=f'{theta_vardict["label"]} {theta_vardict["unit"]}',
        ylabel=f'{var_dict["dtheta"]["label"]} {var_dict["dtheta"]["unit"]}',
    )
    ax1.set_title("Loss function space", loc="left")
    ax1.set_xticks(
        [theta_w, theta_star, theta_fc, n],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$",
            r"$\theta_{\mathrm{fc}}$", r"$n$"],
    )
    ax1.set_yticks(
        [0, -k], [0, r"$\frac{\mathrm{ET}_{\mathrm{max}}}{\Delta z}$"])

    # Update Drydown plot
    ax2.clear()
    # base plot
    ax2.plot(
        t,
        dd_values,
        linewidth=linewidth,
        color=base_color,
        alpha=0.5
    )
    # Drainage plot
    ax2.plot(t_D, dd_values_D, linewidth=linewidth, color=D_color)
    ax2.scatter(t_D[-1], dd_values_D[-1], color=D_color)
    ax2.set_ylim([0.0, n])
    ax2.set_xlim([0, tmax])
    ax2.set(
        xlabel=f'{var_dict["t"]["label"]} {var_dict["t"]["unit"]}',
        ylabel=f'{theta_vardict["label"]} {theta_vardict["unit"]}',
    )
    ax2.set_yticks(
        [theta_w, theta_star, theta_fc, n],
        [r"$\theta_{\mathrm{wp}}$", r"$\theta^{*}$",
            r"$\theta_{\mathrm{fc}}$", r"$n$"],
    )
    ax2.set_title("Drydown space", loc="left")

    fig.tight_layout()
    return ax1, ax2


# Create animation
ani = FuncAnimation(fig, update, frames=200, interval=25, blit=False)
# frames: The total number of frames for the animation.
# interval: The time between frames in milliseconds(e.g., interval=50 for 20 FPS)

# Show animation
plt.show()

# Save animation (optional)
out_dir = "./out"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
# Save the animation as a GIF
ani.save(os.path.join(out_dir, 'output_animation.gif'), writer='pillow', dpi=80)
print(os.path.join(out_dir, 'output_animation.gif'))
