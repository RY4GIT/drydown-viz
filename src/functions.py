# %%
import numpy as np

###################################################################################
# LOSS FUNCTION
###################################################################################


def lossfunc_3stage(theta, ETmax, Ks, beta, n, theta_fc, theta_star, theta_w, q=1.0, z=50):
    """
    Calculate the loss value (-dtheta/dt) corresponding to the soil moisture level (theta).

    Parameters:
    - theta: Current soil moisture level.
    - ETmax: Maximum evapotranspiration rate.
    - Ks: Saturated hydraulic conductivity.
    - beta: Parameter b.
    - n: Porosity.
    - theta_fc: Field capacity saturation.
    - theta_star: Critical soil moisture content.
    - theta_w: Wilting point soil moisture content.
    - q: nonlinearity parameter
    - z: Soil thickness in mm (default is 50 mm).

    Returns:
    - loss: Calculated loss based on theta's range.
    """

    if theta >= theta_fc:
        # Drainage stage
        return loss_drainage(theta=theta, Ks=Ks, beta=beta, n=n, theta_fc=theta_fc, z=z) + loss_ET_stageI(ETmax=ETmax, z=z)

    elif theta_star <= theta < theta_fc:
        # Stage I ET
        return loss_ET_stageI(ETmax=ETmax, z=z)

    elif theta_w <= theta < theta_star:
        # Stage II ET
        return loss_ET_stageII(theta=theta, ETmax=ETmax, theta_star=theta_star, theta_w=theta_w, q=q, z=z)

    elif theta < theta_w and theta >= 0:
        # Below the wilting point, minimal or no loss
        return 0

    else:
        # Handle any negative theta as 0 loss
        return np.nan


def loss_drainage(theta, Ks, beta, theta_fc, n, z=50):
    return -(Ks / z) * (
        (np.exp(beta * (theta - theta_fc)) - 1) /
        (np.exp(beta * (n - theta_fc)) - 1)
    )


def loss_ET_stageI(ETmax, z=50.0):
    return -(ETmax / z)


def loss_ET_stageII(theta, ETmax, theta_star, theta_w, q=1.0, z=50.0):
    return -(ETmax / z) * ((theta - theta_w) / (theta_star - theta_w)) ** q


###################################################################################
# DRYDOWN MODELS
###################################################################################


def theta_3stage(
    t,
    ETmax,
    theta_0,
    theta_w,
    theta_star,
    theta_fc,
    Ks,
    beta,
    n,
    q=1.0,
    z=50.0,
    penalty=0.5,
):
    """
    Calculate the drydown curve for soil moisture over time using non-linear plant stress model.

    Parameters:
        t (int): Timestep, in day.
        z (float): Soil thicness in mm. Default is 50 mm
        ETmax (float): Maximum evapotranpisration rate in mm/day.
        q (float): Degree of non-linearity in the soil moisture response.
        theta_0 (float): The initial soil moisture after precipitation, in m3/m3
        theta_star (float, optional): Critical soil moisture content, equal to s_star * porosity, in m3/m3
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3

    Returns:
        float: Rate of change in soil moisture (dtheta/dt) for the given timestep, in m3/m3/day.
    """
    # ______________________________________________________________
    if (theta_0 <= theta_w) or (theta_star <= theta_w) or (theta_star >= theta_fc):
        return np.full_like(t, penalty)

    if q <= 0:
        return np.full_like(t, penalty)

    # ______________________________________________________________
    # Get the time range

    # If theta_0 is in drainage stage (larger than theta_fc),
    # find the time to reach theta_fc and theta_star
    if theta_0 >= theta_fc:
        try:
            t_fc = find_t_fc(
                ETmax=ETmax, Ks=Ks, beta=beta, theta_fc=theta_fc, theta_0=theta_0, n=n, z=z)
            t_star = find_t_star(ETmax=ETmax, theta_0=theta_fc,
                                 theta_star=theta_star) + t_fc
        except Exception as e:
            # Catch any unexpected errors and set t_fc to NaN
            print(f"Error in finding t_fc: {e}")
            return np.full_like(t, penalty)

    # If theta_0 is in stage I (smaller than theta_fc but larger than theta_star),
    # find the time to reach theta_star
    elif (theta_0 > theta_star) & (theta_0 < theta_fc):
        try:
            t_fc = 0
            t_star = find_t_star(
                ETmax=ETmax, theta_0=theta_0, theta_star=theta_star)
        except Exception as e:
            # Catch any unexpected errors and set t_fc to NaN
            print(f"Error in finding t_fc: {e}")
            return np.full_like(t, penalty)

    # If theta_0 is smaller than theta_star, theta_0 is in stage II ET at t=0
    else:
        t_fc = 0
        t_star = 0

    # ______________________________________________________________
    # Get the theta_0 for each stage
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
                t=t, ETmax=ETmax, Ks=Ks, beta=beta, theta_fc=theta_fc, theta_0=theta_0, n=n),
            lambda t: theta_ET_stageI(
                t=t, ETmax=ETmax, theta_0=theta_0_i, t_fc=t_fc),
            lambda t: theta_ET_expmodel(
                t=t,
                ETmax=ETmax,
                theta_0=theta_0_ii,
                theta_star=theta_star,
                theta_w=theta_w,
                t_star=t_star,
            ),
        ],
    )

    return theta_t


def find_t_fc(ETmax, Ks, beta, theta_fc, theta_0, n, z=50):
    eta = ETmax / z
    m = Ks / (z * (np.exp(beta * (n - theta_fc)) - 1))

    # Second equation for t_sfc
    ln_term = np.log((eta - m + m * np.exp(beta * (theta_0 - theta_fc))) / eta)
    t_fc = (1 / (beta * (m - eta))) * (beta * (theta_fc - theta_0) + ln_term)

    return t_fc


def find_t_star(ETmax, theta_0, theta_star, z=50):
    k = (
        ETmax / z
    )  # Constant term. Convert ETmax to maximum dtheta/dt rate from a unit volume of soil

    # Time it takes from theta_0 to theta_star
    t_star = (theta_0 - theta_star) / k
    return t_star


def theta_drainage(t, ETmax, Ks, beta, theta_fc, theta_0, n, z=50):
    """
        Compute the value of s(t) based on the analytical solution of the differential equation.

        Parameters:
        - t: Time variable (scalar or NumPy array)
        - ETmax: Potential Evapotranspiration (scalar)
        - def compute_s(t, ETmax, Ks, b, theta_fc, theta_0):
    : Hydraulic conductivity constant (scalar)
        - beta: Parameter b in the equation (scalar)
        - theta_fc: Field capacity saturation (scalar)
        - theta_0: Initial saturation at time t=0 (scalar)
        - n: porosity

        Returns:
        - theta_t: Value(s) of s at time t (scalar or NumPy array matching the shape of t)
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


def theta_ET_stageI(t, ETmax, theta_0, t_fc=0.0, z=50):
    return theta_0 - (ETmax / z) * (t - t_fc)


def theta_ET_qmodel(t, q, ETmax, theta_0, theta_star, theta_w, z=50.0, t_star=0.0):

    k = (
        ETmax / z
    )  # Constant term. Convert ETmax to maximum dtheta/dt rate from a unit volume of soil

    b = (theta_0 - theta_w) ** (1 - q)

    a = (1 - q) / ((theta_star - theta_w) ** q)

    base = -k * a * (t - t_star) + b

    return base ** (1 / (1 - q)) + theta_w


def theta_ET_expmodel(t, ETmax, theta_0, theta_star, theta_w, z=50.0, t_star=0.0):
    """
    Calculate the drydown curve for soil moisture over time using non-linear plant stress model.

    Parameters:
        t (int): Timestep, in day.
        z (float): Soil thicness in mm. Default is 50 mm
        ETmax (float): Maximum evapotranpisration rate in mm/day.
        theta_0 (float): The initial soil moisture after precipitation, in m3/m3
        theta_star (float, optional): Critical soil moisture content, equal to s_star * porosity, in m3/m3
        theta_w (float, optional): Wilting point soil moisture content, equal to s_star * porosity, in m3/m3

    Returns:
        float: Rate of change in soil moisture (dtheta/dt) for the given timestep, in m3/m3/day.
    """

    tau = z * (theta_star - theta_w) / ETmax

    if theta_0 > theta_star:
        theta_0_ii = theta_star
    else:
        theta_0_ii = theta_0

    return (theta_0_ii - theta_w) * np.exp(-(t - t_star) / tau) + theta_w
