# %%
import os
import matplotlib.pyplot as plt
import numpy as np

from functions import (
    theta_3stage,
    lossfunc_3stage,
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

Ks = 12
beta = 12
n = 0.7
theta_fc = 0.50
theta_star = 0.45
theta_w = 0.05

# Define variables
theta = np.arange(0, n, 1e-03)
theta_0 = n

tmax = 15
t_star = (theta_fc - theta_star) / k
t_before_star = np.arange(0, t_star, 1e-03)
t = np.arange(0, tmax, 1e-03)

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


# %% Plot
fig = plt.figure(figsize=(8, 4))
linewidth = 3

# _________________________________________________________________________
color = "#5ab4ac"

# _________________________________________________________________________
# Loss function plot
ax1 = fig.add_subplot(1, 2, 1)

loss_values = [lossfunc_3stage(
    theta=theta_t, q=q, ETmax=ETmax, Ks=Ks, beta=beta, n=n, theta_fc=theta_fc, theta_star=theta_star, theta_w=theta_w, z=z) for theta_t in theta]

# base plot
ax1.plot(
    theta,
    loss_values,
    linewidth=linewidth,
    color=color,
)

# TODO: Drainage phase
# ax1.plot(
#     theta,
#     loss_values,
#     linewidth=linewidth,
#     color=color,
# )

ax1.invert_yaxis()
ax1.set_xlim([0.0, n])
ax1.set(
    xlabel=f'{theta_vardict["label"]} {theta_vardict["unit"]}',
    ylabel=f'{dtheta_vardict["label"]} {dtheta_vardict["unit"]}',
)
ax1.set_title(
    label="Loss function space", loc="left"
)
ax1.set_xticks(
    [theta_w, theta_star, theta_fc, n],
    [r"$\theta_{\mathrm{wp}}$", r"$\theta_{*}$",
        r"$\theta_{\mathrm{fc}}$", r"$n$"],
)
ax1.set_yticks([0, -k], [0, r"$\frac{\mathrm{ET}_{\mathrm{max}}}{\Delta z}$"])

# _________________________________________________________________________
# Drydown plot

ax2 = fig.add_subplot(1, 2, 2)

dd_values = [theta_3stage(
    t=t_i, q=q, ETmax=ETmax, theta_0=theta_0, theta_w=theta_w, theta_star=theta_star, theta_fc=theta_fc, Ks=Ks, beta=beta, n=n) for t_i in t]

# base plot
ax2.plot(
    t,
    dd_values,
    linewidth=linewidth,
    color=color,
)

ax2.set_ylim([0.0, n])
ax2.set_xlim([0, tmax])
ax2.set(
    xlabel=f'{var_dict["t"]["label"]} {var_dict["t"]["unit"]}',
    ylabel=f'{theta_vardict["label"]} {theta_vardict["unit"]}',
)
ax2.set_title(label="Drydown space", loc="left")
ax2.set_yticks(
    [theta_w, theta_star, theta_fc, n],
    [r"$\theta_{\mathrm{wp}}$", r"$\theta_*$",
        r"$\theta_{\mathrm{fc}}$", r"$n$"],
)
fig.tight_layout()

# %%
len(dd_values)
len(t)
# %% Output

out_dir = r".\out"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    print(f"Created dir: {out_dir}")


# %%
fig.savefig(os.path.join(out_dir, "test.png"),
            dpi=600, bbox_inches="tight")

# %%
# TODO: Make an animation for drainage and ET phase
