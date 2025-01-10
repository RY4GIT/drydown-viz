import numpy as np

###################################################################################
###################################################################################
# LOSS FUNCTION
###################################################################################
###################################################################################


def lossfunc_3stage(
    theta, ETmax, Ks, beta, n, theta_fc, theta_star, theta_w, q=1.0, z=50
):
    """
    Computes the loss value (-dtheta/dt) for a given soil moisture level (theta)
    based on a three-stage soil moisture loss model.

    This model assumes three distinct stages of soil moisture dynamics:
    1. **Drainage**: Occurs when soil moisture is above the field capacity (theta > theta_fc).
        Loss is dominated by gravitational drainage.
    2. **Stage I ET (Evapotranspiration)**: Occurs when soil moisture is between the field capacity
        and the critical moisture level (theta_fc >= theta > theta_star).
        Loss is dominated by atmospheric demand, with potential evapotranspiration at its maximum (ETmax).
    3. **Stage II ET (Evapotranspiration)**: Occurs when soil moisture is between the critical moisture
        level and the wilting point (theta_star >= theta > theta_w). Loss is limited by soil moisture availability,
        transitioning to a non-linear reduction in evapotranspiration.

    Hygroscopic water (soil moisture below theta_w) is neglected in this model.

    Parameters:
    ----------
    theta : float
        Current soil moisture level (m³/m³).
    ETmax : float
        Maximum evapotranspiration rate (mm/day).
    Ks : float
        Saturated hydraulic conductivity (mm/day).
    beta : float
        Soil moisture stress parameter (dimensionless).
    n : float
        Soil porosity (m³/m³).
    theta_fc : float
        Field capacity, the soil moisture level at which drainage effectively ceases (m³/m³).
    theta_star : float
        Critical soil moisture content for transition between stages (m³/m³).
    theta_w : float
        Wilting point soil moisture content, below which plants cannot extract water (m³/m³).
    q : float, optional
        Nonlinearity parameter (dimensionless), default is 1.0.
    z : float, optional
        Soil thickness (mm), default is 50 mm.

    Returns:
    -------
    float
        The calculated loss value (-dtheta/dt) (m³/m³/day).
    """

    if theta >= theta_fc:
        # Drainage stage
        return loss_drainage(
            theta=theta, Ks=Ks, beta=beta, n=n, theta_fc=theta_fc, z=z
        ) + loss_ET_stageI(ETmax=ETmax, z=z)

    elif theta_star <= theta < theta_fc:
        # Stage I ET
        return loss_ET_stageI(ETmax=ETmax, z=z)

    elif theta_w <= theta < theta_star:
        # Stage II ET
        return loss_ET_stageII(
            theta=theta, ETmax=ETmax, theta_star=theta_star, theta_w=theta_w, q=q, z=z
        )

    elif theta < theta_w and theta >= 0:
        # Below the wilting point, minimal or no loss
        return 0

    else:
        # Handle any negative theta as 0 loss
        return np.nan


###################################################################################
# LOSS MODELS (EACH STAGE)


def loss_drainage(theta, Ks, beta, theta_fc, n, z=50):
    """
    Computes the loss value (-dtheta/dt) for a given soil moisture level (theta) in Drainage stage.
    """
    return -(Ks / z) * (
        (np.exp(beta * (theta - theta_fc)) - 1) / (np.exp(beta * (n - theta_fc)) - 1)
    )


def loss_ET_stageI(ETmax, z=50.0):
    """
    Computes the loss value (-dtheta/dt) for a given soil moisture level (theta) in Stage I ET.
    """
    return -(ETmax / z)


def loss_ET_stageII(theta, ETmax, theta_star, theta_w, q=1.0, z=50.0):
    """
    Computes the loss value (-dtheta/dt) for a given soil moisture level (theta) in Stage II ET.
    """
    return -(ETmax / z) * ((theta - theta_w) / (theta_star - theta_w)) ** q


###################################################################################
###################################################################################
# DRYDOWN MODELS
###################################################################################
###################################################################################


def theta_3stage(
    t,
    theta_0,
    ETmax,
    Ks,
    beta,
    n,
    theta_w,
    theta_star,
    theta_fc,
    q=1.0,
    z=50.0,
    penalty=0.5,
):
    """
    Computes the soil moisture (theta) at a given timestep (t) based on a three-stage
    drydown curve model for soil moisture dynamics.

    This model incorporates three stages of soil moisture dynamics:
    1. **Drainage**: Occurs when soil moisture (theta) is greater than the field capacity (theta > theta_fc).
        Loss is dominated by gravitational drainage.
    2. **Stage I ET (Evapotranspiration)**: Occurs when theta lies between field capacity and the
        critical soil moisture content (theta_fc >= theta > theta_star). Loss is primarily determined
        by atmospheric demand, with evapotranspiration occurring at its maximum rate (ETmax).
    3. **Stage II ET (Evapotranspiration)**: Occurs when theta is between the critical soil moisture
        content and the wilting point (theta_star >= theta > theta_w). Loss is constrained by soil moisture
        availability, with evapotranspiration reducing non-linearly.

    Hygroscopic water (below wilting point, theta_w) is excluded from this model.

    Parameters:
    ----------
    t : int
        Timestep, in days.
    theta_0 : float
        Initial soil moisture after precipitation, in m³/m³.
    ETmax : float
        Maximum evapotranspiration rate, in mm/day.
    Ks : float
        Saturated hydraulic conductivity, in mm/day.
    beta : float
        Soil moisture stress parameter (dimensionless).
    n : float
        Soil porosity, in m³/m³.

    theta_w : float
        Wilting point soil moisture content, equal to s_w * porosity, in m³/m³.
    theta_star : float
        Critical soil moisture content, equal to s_star * porosity, in m³/m³.
    theta_fc : float
        Field capacity soil moisture content, in m³/m³.
    q : float, optional
        Degree of non-linearity in the soil moisture response (dimensionless). Default is 1.0.
    z : float, optional
        Soil thickness, in mm. Default is 50 mm.
    penalty : float, optional
        Penalty parameter for failed estimation. Default is 0.5.

    Returns:
    -------
    float
        Soil moisture content (theta) at the given timestep, in m³/m³.
    """

    # ______________________________________________________________\
    # Parameter checks
    if (theta_0 <= theta_w) or (theta_star <= theta_w) or (theta_star >= theta_fc):
        return np.full_like(t, penalty)

    if q <= 0:
        return np.full_like(t, penalty)

    # ______________________________________________________________
    # Calculate Time Ranges for Each Stage

    try:
        if theta_0 >= theta_fc:  # Stage: Drainage
            t_fc = find_t_fc(
                theta_0=theta_0,
                ETmax=ETmax,
                Ks=Ks,
                beta=beta,
                n=n,
                theta_fc=theta_fc,
                z=z,
            )
            t_star = (
                find_t_star(theta_0=theta_fc, ETmax=ETmax, theta_star=theta_star) + t_fc
            )

        elif theta_0 > theta_star:  # Stage: Stage I ET
            t_fc = 0
            t_star = find_t_star(theta_0=theta_0, ETmax=ETmax, theta_star=theta_star)

        else:  # Stage: Stage II ET
            t_fc = 0
            t_star = 0

    except Exception as e:
        print(f"Error in calculating time ranges: {e}")
        return np.full_like(t, penalty)

    # ______________________________________________________________
    # Get the initial soil mositure values entring into each stage
    if theta_0 > theta_fc:
        theta_0_i = theta_fc
    else:
        theta_0_i = theta_0

    if theta_0 > theta_star:
        theta_0_ii = theta_star
    else:
        theta_0_ii = theta_0

    # ______________________________________________________________
    # Use np.piecewise to apply different functions based on the conditions
    theta_t = np.piecewise(
        t,
        [t < t_fc, (t >= t_fc) & (t < t_star), t >= t_star],
        [
            lambda t: theta_drainage(
                t=t,
                ETmax=ETmax,
                Ks=Ks,
                beta=beta,
                theta_fc=theta_fc,
                theta_0=theta_0,
                n=n,
                z=z,
            ),
            lambda t: theta_ET_stageI(
                t=t, ETmax=ETmax, theta_0=theta_0_i, t_fc=t_fc, z=z
            ),
            lambda t: theta_ET_stageII(
                t=t,
                q=q,
                ETmax=ETmax,
                theta_0=theta_0_ii,
                theta_star=theta_star,
                theta_w=theta_w,
                t_star=t_star,
                z=z,
            ),
        ],
    )

    return theta_t


###################################################################################
# TIME FUNCTIONS


def find_t_fc(theta_0, ETmax, Ks, beta, n, theta_fc, z=50):
    """This function computes the time (t_fc) required for the soil moisture level (theta_0)
    to decrease to the field capacity soil moisture content (theta_fc) during Drainage Stage.
    """
    eta = ETmax / z
    m = Ks / (z * (np.exp(beta * (n - theta_fc)) - 1))
    ln_term = np.log((eta - m + m * np.exp(beta * (theta_0 - theta_fc))) / eta)
    t_fc = (1 / (beta * (m - eta))) * (beta * (theta_fc - theta_0) + ln_term)

    return t_fc


def find_t_star(theta_0, ETmax, theta_star, z=50):
    """This function computes the time (t_star) required for the soil moisture level (theta_0)
    to decrease to the critical soil moisture content (theta_star) during Stage I evapotranspiration.
    """
    k = ETmax / z
    t_star = (theta_0 - theta_star) / k
    return t_star


###################################################################################
# DRYDOWN MODELS (EACH STAGE)


def theta_drainage(t, theta_0, ETmax, Ks, n, beta, theta_fc, z=50):
    """
    Calculate the drydown curve for soil moisture over time in Drainage Stage.
    """

    # Compute constants
    En = np.exp(beta * (n - theta_fc))
    E0 = np.exp(beta * (theta_0 - theta_fc))

    C = En - 1
    q = (beta * Ks) / (C * z)
    p = (beta / z) * (ETmax - (Ks / C))

    # Compute numerator and denominator for E(t)
    numerator = p * E0 * np.exp(-p * t)
    denominator = p + q * E0 * (1 - np.exp(-p * t))

    # Compute s(t)
    fraction = numerator / denominator
    theta_t = theta_fc + (1 / beta) * np.log(fraction)

    return theta_t


def theta_ET_stageI(t, theta_0, ETmax, t_fc=0.0, z=50):
    """
    Calculate the drydown curve for soil moisture over time in Stage II ET.
    """
    return theta_0 - (ETmax / z) * (t - t_fc)


def theta_ET_stageII(t, theta_0, ETmax, theta_star, theta_w, q=1.0, z=50.0, t_star=0.0):
    """
    Calculate the drydown curve for soil moisture over time in Stage II ET.
    """

    if q == 1.0:
        return theta_ET_expmodel(
            t=t,
            theta_0=theta_0,
            ETmax=ETmax,
            theta_star=theta_star,
            theta_w=theta_w,
            z=z,
            t_star=t_star,
        )
    else:
        return theta_ET_qmodel(
            t=t,
            theta_0=theta_0,
            ETmax=ETmax,
            theta_star=theta_star,
            theta_w=theta_w,
            q=q,
            z=z,
            t_star=t_star,
        )


def theta_ET_expmodel(t, theta_0, ETmax, theta_star, theta_w, z=50.0, t_star=0.0):
    """
    Calculate the drydown curve for soil moisture over time in Stage II ET using linear plant stress model (q=1).
    """

    tau = z * (theta_star - theta_w) / ETmax

    return (theta_0 - theta_w) * np.exp(-(t - t_star) / tau) + theta_w


def theta_ET_qmodel(t, theta_0, ETmax, theta_star, theta_w, q, z=50.0, t_star=0.0):
    """
    Calculate the drydown curve for soil moisture over time in Stage II ET using non-linear plant stress model (q!=1).
    """

    k = ETmax / z

    b = (theta_0 - theta_w) ** (1 - q)

    a = (1 - q) / ((theta_star - theta_w) ** q)

    base = -k * a * (t - t_star) + b

    return base ** (1 / (1 - q)) + theta_w
