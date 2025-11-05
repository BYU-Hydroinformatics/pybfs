# -*- coding: utf-8 -*-
"""PyBFS - Python Baseflow Separation

A Python implementation of Baseflow Separation algorithms for hydrological analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import math


def sur_z(lb, a, ws, por, ss):
    """Calculates saturated thickness of the surface reservoir (zs)

    Uses the quadratic formula to solve for saturated thickness based on
    surface storage and basin geometry. If the discriminant is negative,
    returns the maximum possible saturated thickness.

    Parameters
    ----------
    lb : float
        Basin length (m)
    a : float
        Alpha parameter - shape parameter controlling reservoir geometry
    ws : float
        Surface width (m)
    por : float
        Porosity (dimensionless, 0-1)
    ss : float
        Surface storage (m³)

    Returns
    -------
    float
        Saturated thickness (m)

    Notes
    -----
    The function solves: ss = lb * por * (2 * ws * zs - zs²/a)
    """
    a1 = 1 / (2*a)
    b1 = -2 * ws
    c1 = ss / (lb * por)

    discriminant = b1 ** 2 - 4 * a1 * c1

    if discriminant < 0:
        return ws * a
    else:
        return (-b1 - math.sqrt(discriminant)) / (2 * a1)


def sur_store(lb, a, ws, por, zs):
    """Calculates surface storage (ss) in the surface reservoir

    Computes the volume of water stored in the surface reservoir based on
    basin geometry, porosity, and saturated thickness.

    Parameters
    ----------
    lb : float
        Basin length (m)
    a : float
        Alpha parameter - shape parameter controlling reservoir geometry
    ws : float
        Surface width (m)
    por : float
        Porosity (dimensionless, 0-1)
    zs : float
        Saturated thickness (m)

    Returns
    -------
    float
        Surface storage volume (m³)

    Notes
    -----
    Formula: ss = lb * por * (2 * ws * zs - zs²/a)
    """
    z = min(ws * a, zs)
    result = lb * (2 * ws * zs - zs**2 / a) * por
    return result


def sur_q(lb, a, ks, z):
    """Calculates the surface discharge from reservoir (Qs)

    Computes the outflow rate from the surface reservoir based on Darcy's law
    and basin geometry.

    Parameters
    ----------
    lb : float
        Basin length (m)
    a : float
        Alpha parameter - shape parameter controlling reservoir geometry
    ks : float
        Surface hydraulic conductivity (m/day)
    z : float
        Water surface elevation (m)

    Returns
    -------
    float
        Surface discharge rate (m³/day)

    Notes
    -----
    Formula: Qs = 2 * lb * z * a * ks
    """
    return 2 * lb * z * a * ks


def dir_q(lb, a, z, i):
    """Calculates the direct runoff (Qd) from the surface reservoir

    Computes the direct runoff that occurs when precipitation falls on
    saturated areas and immediately becomes streamflow without infiltration.

    Parameters
    ----------
    lb : float
        Basin length (m)
    a : float
        Alpha parameter - shape parameter controlling reservoir geometry
    z : float
        Water surface elevation (m)
    i : float
        Impulse/precipitation intensity (m/day)

    Returns
    -------
    float
        Direct runoff rate (m³/day)

    Notes
    -----
    Formula: Qd = 2 * lb * z / a * i
    Represents runoff from saturated contributing areas
    """
    return 2 * lb * z / a * i


def infiltration(lb, ws, ks, a, zs, i):
    """Calculates the infiltration rate (I) from surface to subsurface

    Computes the rate at which water infiltrates from the surface reservoir
    into the subsurface. Infiltration is limited by either the precipitation
    rate or the hydraulic conductivity.

    Parameters
    ----------
    lb : float
        Basin length (m)
    ws : float
        Surface width (m)
    ks : float
        Surface hydraulic conductivity (m/day)
    a : float
        Alpha parameter - shape parameter controlling reservoir geometry
    zs : float
        Saturated thickness (m)
    i : float
        Impulse/precipitation intensity (m/day)

    Returns
    -------
    float
        Infiltration rate (m³/day)

    Notes
    -----
    Formula: I = 2 * lb * (ws - zs/a) * min(i, ks)
    Infiltration occurs in unsaturated areas and is rate-limited
    """
    return 2 * lb * (ws - zs / a) * min(i, ks)


def recharge(lb, xb, ws, kz, zs, por):
    """Calculates recharge rate (R) from surface to base reservoir

    Computes the vertical recharge rate from the surface reservoir to the
    base (groundwater) reservoir. Recharge is limited by either the available
    water or the vertical hydraulic conductivity.

    Parameters
    ----------
    lb : float
        Basin length (m)
    xb : float
        Longitudinal location of base water level intersection (m)
    ws : float
        Surface width (m)
    kz : float
        Vertical hydraulic conductivity (m/day)
    zs : float
        Saturated thickness (m)
    por : float
        Porosity (dimensionless, 0-1)

    Returns
    -------
    float
        Recharge rate (m³/day)

    Notes
    -----
    Formula: R = (lb - xb) * 2 * ws * min(zs * por, kz)
    Represents water moving from surface to groundwater storage
    """
    return (lb - xb) * 2 * ws * min(zs * por, kz)


def get_values_for_site(df, site_no):
    """Extracts site-specific parameters from a DataFrame

    Searches a parameter DataFrame and extracts basin characteristics,
    groundwater hydraulic parameters, and flow metrics for a specified site.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing site parameters with columns site_no (site identification
        number), basin characteristics (AREA, Lb, X1, Wb, POR), groundwater hydraulics
        (ALPHA, BETA, Ks, Kb, Kz), and flow metrics (Qthresh, Rs, Rb1, Rb2, Prec, Frac4Rise)
    site_no : int
        Site identification number to search for

    Returns
    -------
    tuple of (basin_char, gw_hyd, flow)
        basin_char : list of [AREA, Lb, X1, Wb, POR]
        gw_hyd : list of [ALPHA, BETA, Ks, Kb, Kz]
        flow : list of [Qthresh, Rs, Rb1, Rb2, Prec, Frac4Rise]

    Examples
    --------
    >>> params_df = pd.read_csv('site_parameters.csv')
    >>> basin, gw, flow = get_values_for_site(params_df, 2312200)
    """
    # Define the column groups
    basin_char_columns = ["AREA", "Lb", "X1", "Wb", "POR"]
    gw_hyd_columns = ["ALPHA", "BETA", "Ks", "Kb", "Kz"]
    flow_columns = ["Qthresh", "Rs", "Rb1", "Rb2", "Prec", "Frac4Rise"]

    # Initialize the dictionaries
    basin_char_dict = {}
    gw_hyd_dict = {}
    flow_dict = {}

    # Iterate over the rows of the dataframe
    for idx, row in df.iterrows():
        site_no_key = row['site_no']

        # Assign the respective columns to the dictionaries
        basin_char_dict[site_no_key] = list(row[basin_char_columns].values)
        gw_hyd_dict[site_no_key] = list(row[gw_hyd_columns].values)
        flow_dict[site_no_key] = list(row[flow_columns].values)

    # Get the values for the specific site number
    basin_char = basin_char_dict.get(site_no)
    gw_hyd = gw_hyd_dict.get(site_no)
    flow = flow_dict.get(site_no)

    return basin_char, gw_hyd, flow


def base_table(lb, x1, wb, b, kb, q, por):
    """Generates baseflow table relating discharge to storage and water levels

    Creates a lookup table that relates baseflow discharge (Q) to the
    longitudinal location of base water level intersection (Xb), water
    surface elevation (Z), and storage (S) in the base reservoir.

    Parameters
    ----------
    lb : float
        Basin length (m)
    x1 : float
        Initial longitudinal position parameter (m)
    wb : float
        Basin width (m)
    b : float
        Beta parameter, controls the relationship between discharge and water levels
    kb : float
        Base hydraulic conductivity (m/day)
    q : pd.DataFrame
        Streamflow data containing 'Streamflow' column (m³/day)
    por : float
        Porosity (dimensionless, 0-1)

    Returns
    -------
    pd.DataFrame
        Baseflow table with columns Xb (longitudinal location m), Z (water surface
        elevation m), S (storage m³), and Q (discharge m³/day)

    Notes
    -----
    The table is generated using log-spaced discharge values ranging from
    minimum to maximum observed streamflow. Used for lookup during PyBFS simulation.
    """
    qin = np.array(q['Streamflow'])
    # Remove NaNs and extract positive flow values
    tmp_q = qin[~np.isnan(qin)]
    tmp_q = tmp_q[tmp_q > 0]

    # Define the log range of flow values
    tmp_range = np.log10([np.min(tmp_q) / 10, np.max(tmp_q)])

    # Create log-spaced discharge values
    qq = np.logspace(tmp_range[0], tmp_range[1], num=1000)
    qq[0] = 0

    # Calculate z based on the value of b
    if b != 0.5:
        z = ((qq * x1 * (2 * b - 1)) / (wb * kb * b ** 2)) ** (b / (2 * b - 1))
    else:
        z = np.exp(qq * 2 * x1 / (wb * kb))

    # Compute x and storage s
    x = x1 * z ** (1 / b)
    s = wb * por * ((1 / x1 ** b) * (1 / (b + 1)) * x ** (b + 1) + (lb - x) * z)

    # Construct the table and round values
    BT = pd.DataFrame({
        'Xb': np.round(x, 5),
        'Z':  np.round(z, 5),
        'S':  np.round(s, 5),
        'Q':  np.round(qq, 5)
    })

    # Filter out rows where x < lb
    BT = BT[BT['Xb'] < lb]

    return BT


def PyBFS(streamflow, SBT, basin_char, gw_hyd, flow):
    """Main PyBFS function for baseflow separation

    Performs physically-based baseflow separation using a coupled surface-subsurface
    reservoir model. Separates total streamflow into three components: surface flow,
    baseflow, and direct runoff. Uses an iterative approach to estimate precipitation
    impulses needed to match observed streamflow.

    Parameters
    ----------
    streamflow : pd.DataFrame
        DataFrame with columns 'Date' (datetime) and 'Streamflow' (m³/day,
        observed streamflow)
    SBT : pd.DataFrame
        Baseflow table with columns ['Xb','Z','S','Q']. Generated by
        base_table() function
    basin_char : list
        Basin characteristics [area, lb, x1, wb, por] where:

        - area: Basin area (m²)
        - lb: Basin length (m)
        - x1: Initial longitudinal position (m)
        - wb: Basin width (m)
        - por: Porosity (0-1)
    gw_hyd : list
        Groundwater hydraulic parameters [alpha, beta, ks, kb, kz] where:

        - alpha: Surface reservoir shape parameter
        - beta: Base reservoir shape parameter
        - ks: Surface hydraulic conductivity (m/day)
        - kb: Base hydraulic conductivity (m/day)
        - kz: Vertical hydraulic conductivity (m/day)
    flow : list
        Flow metrics [qthresh, rs, rb1, rb2, prec, fr4rise] where:

        - qthresh: Flow threshold
        - rs: Recession slope parameter
        - rb1, rb2: Baseflow recession parameters
        - prec: Precision threshold
        - fr4rise: Fraction for rise detection

    Returns
    -------
    pd.DataFrame
        Results DataFrame with columns:

        - Date: Date of observation
        - Qob: Observed streamflow (m³/day)
        - Qsim: Simulated total streamflow (m³/day)
        - SurfaceFlow: Surface flow component (m³/day)
        - Baseflow: Baseflow component (m³/day)
        - DirectRunoff: Direct runoff component (m³/day)
        - X: Longitudinal location of base water level (m)
        - Eta: Streamflow residual (m³/day)
        - StSur: Surface storage (m³)
        - StBase: Base storage (m³)
        - Impulse.L: Estimated precipitation (m/day)
        - Zs.L: Surface water elevation (m)
        - Zb.L: Base water elevation (m)
        - Infil: Infiltration rate (m³/day)
        - Rech: Recharge rate (m³/day)
        - RecessCount.T: Recession day counter
        - AdjPctEr: Adjusted percent error
        - Weight: Error weighting factor

    Notes
    -----
    The model uses a two-reservoir system where surface and base reservoirs
    interact through infiltration and recharge. The algorithm iteratively
    estimates precipitation impulses during rising limbs to minimize streamflow
    residuals (Eta).

    Examples
    --------
    >>> streamflow_data = pd.read_csv('streamflow.csv')
    >>> params_df = pd.read_csv('site_parameters.csv')
    >>> basin_char, gw_hyd, flow = get_values_for_site(params_df, site_no)
    >>> SBT = base_table(lb, x1, wb, beta, kb, streamflow_data, por)
    >>> results = PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
    """
    date = pd.to_datetime(streamflow["Date"])
    qin = np.array(streamflow['Streamflow'])
    qmean = np.nanmean(qin)

    timestep = 'day'
    error_basis = 'total'

    # Error tolerance used sequentially to refine impulse
    ifact = [2, 1.1]
    # Number of time steps
    p = len(qin)

    # calculates the change in streamflow (dq) between consecutive time steps. This change helps identify recessional periods (when streamflow is decreasing).
    dq = np.zeros(p)
    if timestep == 'day':
        dq[1:] = qin[1:] - qin[:-1]
    elif timestep == 'hour':
        for y in range(24, p):
            dq[y] = qin[y] - np.nanmax(qin[(y-24):y])

    #basin characteristics
    area, lb, x1, wb, por, ws = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4], basin_char[3] / 2

    # Groundwater hydraulic parameters
    alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]

    # Flow metrics
    qthresh, rs, rb1, rb2, prec, fr4rise = flow[0], flow[1], flow[2], flow[3], flow[4], flow[5]

    #dqfr represents the fractional change in streamflow. #It's used to determine the nature of the change in streamflow relative to the current flow.
    dqfr = dq / qin
    dqfr[(dq == 0) & (qin == 0)] = 0
    dqfr[(dq < 0) & (qin == 0)] = 1

    rise = (dqfr > fr4rise) & (dq > prec)
    rise[np.isnan(rise)] = False

    recess = (dqfr <= fr4rise) | (dq < prec)
    recess[np.isnan(recess)] = False

    #recess_day: An array counting the number of consecutive recession periods up to each time step
    recess_day = np.cumsum(recess) - np.maximum.accumulate((~recess).astype(int) * np.cumsum(recess))

    # Output variables
    X = np.full(p, np.nan)  #LONGITUDINAL LOCATION OF BASE WATER LEVEL INTERSECTION WITH SURFACE, xb
    qcomp = np.full((p, 3), np.nan) ##THREE FLOW COMPONENTS: surface flow; base flow; direct runoff from saturated areas
    ETA = np.full(p, np.nan)  #STATE DISTURBANCES (POSITIVE VALUES REPRESENT INPUTS)
    I = np.full(p, np.nan) #PRECIPITATION CALCULATED FROM eta
    Z = np.full((p, 2), np.nan) #WATER SURFACE ELEVATION OF SURFACE (CHANNEL IS DATUM) AND BASE (BASIN OUTLET IS DATUM), ZS and Zb
    ST = np.full((p, 2), np.nan) #STORAGE, surface and base
    EXC = np.full((p, 2), np.nan) #EXCHANGES, INFILTRATION AND RECHARGE


    #CHECK PARAMETERS, END PROCESS IF PARAMETERS ARE BAD
    if np.any(np.array([lb, x1, wb, alpha, beta, ks, kb, ks, por, qthresh, -rs, -rb1, -rb2, prec, fr4rise]) < 0):
        ts = 10 * p # 10 * p is considered invalid
        print('Negative parameter(s)')
    if lb * wb > area:
        ts = 10 * p
        print('lb x wb > area')
    if np.any(np.isnan(SBT)):
        ts = 10 * p
        print('Cannot calculate discharge for base parameters')

    #basin characteristics
    area, lb, x1, wb, por, ws = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4], basin_char[3] / 2

    # Groundwater hydraulic parameters
    alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]

    # Flow metrics
    qthresh, rs, rb1, rb2, prec, fr4rise = flow[0], flow[1], flow[2], flow[3], flow[4], flow[5]

    ts = 0  #INITIAL TIME STEP
    stts = ts  #STARTING TIME STEP, stts, FOR ERROR CALCULATION


    #This code iterates through the qin array starting from a time step ts, skipping over any NaN values.
    #When it finds the first valid (non-NaN) value, it stores that time step in the variable sttts.
    while np.isnan(qin[ts]):
        ts += 1
        stts = ts #
    # Initialize Variables
    ts_ini = True
    # Initialize variables
    qb_in = min(qin[ts], qmean)
    qb_en = np.nan
    idx = (SBT["Q"] <= qb_in).sum()

    # Adjust for zero-based indexing and avoid out-of-bounds errors
    xb_in = SBT["Xb"].iloc[idx - 1] if idx > 0 else np.nan
    zb_in = SBT["Z"].iloc[idx - 1] if idx > 0 else np.nan
    sb_in = SBT["S"].iloc[idx - 1] if idx > 0 else np.nan

    # Calculate available base storage
    sba = np.max(SBT['S']) - sb_in

    # Surface flow and other calculations
    qs_in = max(0, qin[ts] - qb_in)  # Ensure surface flow is non-negative
    zs_in = min(qs_in / (2 * lb * ks * alpha), ws * alpha)  # Saturated thickness of surface reservoir
    ss_in = sur_store(lb, alpha, ws, por, zs_in)  # Surface storage
    ssa = sur_store(lb, alpha, ws, por, ws * alpha) - ss_in  # AVAILABLE SURFACE STORAGE

    # Infiltration and recharge
    infil_in = 0
    rech_in = recharge(lb, xb_in, ws, kz, zs_in, por)

    while ts < p:

        # Initialize Time Step Using State Variables for Previous Time Step if Available
        if not ts_ini:
            xb_in = X[ts - 1]
            zb_in = Z[ts - 1, 1]
            sb_in = ST[ts - 1, 1]

            idx = (SBT["Xb"] <= xb_in).sum()
            qb_in = SBT["Q"].iloc[idx - 1] if idx > 0 else np.nan

            zs_in = Z[ts - 1, 0]
            ss_in = ST[ts - 1, 0]
            qs_in = sur_q(lb, alpha, ks, zs_in)

            # Storage Capacity Available
            ssa = sur_store(lb, alpha, ws, por, ws * alpha) - ss_in
            sba = max(SBT['S']) - sb_in
            rech_in = min(recharge(lb, xb_in, ws, kz, zs_in, por), sba + qb_in)
            qd = 0
            infil_in = 0

            # Impulse (PPT) Needed to Generate Observed Streamflow
            I[ts] = 0  # Set Impulse to Zero
            etaest = max(0, qin[ts] - qb_in - qs_in)

            # Initial Estimate of Impulse Required for Additional Surface Flow
            if ts > 1 and etaest > 0:
                if rise[ts] or rise[ts - 1]:
                    I[ts] = etaest / (2 * lb * ws)
                    zs = zs_in
                    qs = qs_in

                    # Loop to Calculate Additional Impulse Needed to Reduce ETA
                    # Use Progressively Smaller Incremental Changes in Impulse (ifact) for Iterations
                    for x in ifact:
                        i = I[ts]
                        eta = etaest
                        while eta > max(prec, qin[ts] / 100) and i > 0:
                            I[ts] = i
                            etaest = eta
                            i = x * i
                            infil = min(infiltration(lb, ws, ks, alpha, (zs_in + zs) / 2, i), ssa)  # Limit Infiltration to Available Storage
                            ss = max(ss_in + infil - rech_in - qs, 0)  # Update Surface Storage
                            zs = sur_z(lb, alpha, ws, por, ss)
                            qs = sur_q(lb, alpha, ks, zs)
                            qd = dir_q(lb, alpha, zs_in, i) + dir_q(lb, alpha, (zs - zs_in), i / 2) + max(2 * lb * (ws - zs_in / alpha) * (I[ts] - ks), 0)
                            eta = qin[ts] - qs - qd - qb_in

            infil_in = min(infiltration(lb, ws, ks, alpha, zs_in, I[ts]), ssa)  # Close Initial Calculations When Streamflow Record is Available (Not Projection)

        # End of Time Step Calculations
        ss_en = max(ss_in + infil_in - rech_in - qs_in, 0)
        zs_en = sur_z(lb, alpha, ws, por, ss_en)
        qs_en = sur_q(lb, alpha, ks, zs_en)
        infil_en = min(infiltration(lb, ws, ks, alpha, zs_en, I[ts]), ssa)
        rech_en = min(recharge(lb, xb_in, ws, kz, zs_en, por), sba + qb_in)
        sb_en = max(sb_in + rech_en - qb_in, 0)
        idx = max((SBT["S"] < sb_en).sum(), 1) - 1

        # Safely extract the values from the DataFrame
        xb_en = SBT["Xb"].iloc[idx] if 0 <= idx < len(SBT) else np.nan
        zb_en = SBT["Z"].iloc[idx] if 0 <= idx < len(SBT) else np.nan
        qb_en = SBT["Q"].iloc[idx] if 0 <= idx < len(SBT) else np.nan

        # Final Calculations for Time Step
        qcomp[ts, 0] = (qs_in + qs_en) / 2  # Surface Flow
        qcomp[ts, 1] = (qb_in + qb_en) / 2  # Base Flow

        EXC[ts, 0] = (infil_in + infil_en) / 2
        EXC[ts, 1] = (rech_in + rech_en) / 2

        # For Initial Time Step
        if ts_ini:
            ST[ts, :] = [ss_en, sb_en]
            Z[ts, :] = [zs_en, zb_en]

        # For Time Steps When States Are Available for Previous Time Step
        if not ts_ini:
            ST[ts, 0] = max(ST[ts - 1, 0] + EXC[ts, 0] - qcomp[ts, 0] - EXC[ts, 1], 0)
            ST[ts, 0] = min(ST[ts, 0], sur_store(lb, alpha, ws, por, ws * alpha))
            Z[ts, 0] = sur_z(lb, alpha, ws, por, ST[ts, 0])
            ST[ts, 1] = max(ST[ts - 1, 1] + EXC[ts, 1] - qcomp[ts, 1], 0)
            ST[ts, 1] = min(ST[ts, 1], max(SBT['S']))

            idx = max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1
            Z[ts, 1] = SBT['Z'].iloc[idx] if 0 <= idx < len(SBT) else np.nan

            # Direct Runoff includes additional saturated area x half of rainfall (excess after infiltration), and any precipitation that exceeds infiltration rate
            qcomp[ts, 2] = dir_q(lb, alpha, zs_in, I[ts]) + dir_q(lb, alpha, (Z[ts, 0] - zs_in), I[ts] / 2) + max(2 * lb * (ws - zs_in / alpha) * (I[ts] - ks), 0)

        ETA[ts] = qin[ts] - np.sum(qcomp[ts, 0:3])  # Streamflow Residual
        idx = max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1
        X[ts] = SBT['Xb'].iloc[idx] if 0 <= idx < len(SBT) else np.nan

        ts += 1
        ts_ini = False
        #CLOSE CONDITION ts<p

    if error_basis == 'base':
        q4er = qcomp[:, 1]
    elif error_basis == 'total':
        q4er = np.sum(qcomp, axis=1)

    # ADJUSTED PERCENT ERROR
    APE = (q4er - qin) / (qin + prec)
    APE[(qin == 0) & (q4er == 0)] = 0

    # WEIGHT VARIES FROM 0 TO 1 WITH INCREASING LENGTH OF RECESSION
    Weight = 1 - np.exp(rs * recess_day)

    # WEIGHT OF 1 IS ASSIGNED TO OVER PREDICTION
    Weight[APE > 0] = 1

    tmp = pd.DataFrame({'Date': date, 'Qob': qin, 'Qsim': np.sum(qcomp, axis=1), 'SurfaceFlow': qcomp[:, 0], 'Baseflow': qcomp[:, 1], 'DirectRunoff': qcomp[:, 2], 'X': X,'Eta': ETA, 'StSur': ST[:, 0], 'StBase': ST[:, 1], 'Impulse.L': I, 'Zs.L': Z[:, 0], 'Zb.L': Z[:, 1], 'Infil': EXC[:, 0], 'Rech': EXC[:, 1], 'RecessCount.T': recess_day, 'AdjPctEr': APE, 'Weight': Weight})
    tmp = tmp[['Date', 'Qob', 'Qsim', 'SurfaceFlow', 'Baseflow', 'DirectRunoff',  'X','Eta', 'StSur', 'StBase', 'Impulse.L', 'Zs.L', 'Zb.L', 'Infil', 'Rech', 'RecessCount.T', 'AdjPctEr', 'Weight']]
    return tmp


def forecast(streamflow, SBT, basin_char, gw_hyd, flow, initial):
    """Forecast baseflow using PyBFS with provided initial conditions

    Projects baseflow and reservoir states forward in time using the PyBFS model
    without observed streamflow data. Uses initial conditions from a previous
    PyBFS run to initialize the forecast period.

    Parameters
    ----------
    streamflow : pd.DataFrame
        DataFrame with columns:
        - 'date': Datetime for forecast period
        - 'streamflow': Set to NaN (no observations during forecast)
    SBT : pd.DataFrame
        Baseflow table with columns ['Xb','Z','S','Q']
        Same table used in calibration period
    basin_char : list
        Basin characteristics [area, lb, x1, wb, por]
        Same as used in PyBFS()
    gw_hyd : list
        Groundwater hydraulic parameters [alpha, beta, ks, kb, kz]
        Same as used in PyBFS()
    flow : list
        Flow metrics [qthresh, rs, rb1, rb2, prec, fr4rise]
        Same as used in PyBFS()
    initial : tuple
        Initial states from last time step of calibration:
        (Xi, Zbi, Zsi, StBi, StSi, Surflow, Baseflow, Rech)
        - Xi: Longitudinal location of base water level (m)
        - Zbi: Base water elevation (m)
        - Zsi: Surface water elevation (m)
        - StBi: Base storage (m³)
        - StSi: Surface storage (m³)
        - Surflow: Surface flow (m³/day)
        - Baseflow: Baseflow (m³/day)
        - Rech: Recharge rate (m³/day)

    Returns
    -------
    pd.DataFrame
        Forecasted flow components and states with columns:
        - 'Date': Date of forecast
        - 'Baseflow': Forecasted baseflow (m³/day)
        - 'StSur': Surface storage (m³)
        - 'StBase': Base storage (m³)
        - 'Zs': Surface water elevation (m)
        - 'Zb': Base water elevation (m)
        - 'Rech': Recharge rate (m³/day)

    Notes
    -----
    This function assumes no precipitation during the forecast period (Impulse = 0).
    It projects how baseflow and storage will evolve based solely on drainage and
    reservoir interactions. Best used for short-term (days to weeks) forecasts.

    Examples
    --------
    >>> # Run calibration period
    >>> result = PyBFS(streamflow_cal, SBT, basin_char, gw_hyd, flow)
    >>> # Extract initial conditions from last time step
    >>> ini = result.iloc[-1][['X', 'Zb.L', 'Zs.L', 'StBase', 'StSur',
    ...                         'SurfaceFlow', 'Baseflow', 'Rech']]
    >>> # Create forecast DataFrame
    >>> forecast_dates = pd.date_range(start='2018-10-01', end='2018-11-30')
    >>> forecast_df = pd.DataFrame({'date': forecast_dates, 'streamflow': np.nan})
    >>> # Run forecast
    >>> forecast_result = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
    """
    #basin characteristics
    area, lb, x1, wb, por, ws = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4], basin_char[3] / 2

    # Groundwater hydraulic parameters
    alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]

    # Flow metrics
    qthresh, rs, rb1, rb2, prec, fr4rise = flow[0], flow[1], flow[2], flow[3], flow[4], flow[5]

    date = pd.to_datetime(streamflow["date"])
    qin = np.full(len(np.array(streamflow['streamflow'])), np.nan)
    p = len(qin)
    # Output variables
    X = np.full(p, np.nan)  #LONGITUDINAL LOCATION OF BASE WATER LEVEL INTERSECTION WITH SURFACE, xb
    qcomp = np.full((p, 3), np.nan) ##THREE FLOW COMPONENTS: surface flow; base flow; direct runoff from saturated areas
    ETA = np.full(p, np.nan)  ##STATE DISTURBANCES (POSITIVE VALUES REPRESENT INPUTS) [L3]
    I = np.full(p, np.nan) #PRECIPITATION CALCULATED FROM eta
    Z = np.full((p, 2), np.nan) #WATER SURFACE ELEVATION OF SURFACE (CHANNEL IS DATUM) AND BASE (BASIN OUTLET IS DATUM), ZS and Zb
    ST = np.full((p, 2), np.nan) #STORAGE, surface and base
    EXC = np.full((p, 2), np.nan) #EXCHANGES, INFILTRATION AND RECHARGE# Initialize Variables

    Xi, Zbi, Zsi, StBi, StSi, Surflow, Baseflow, Rech = initial

    ts=0
    ts_ini = True

    while ts < p:
        if ts_ini:
            xb_in = Xi
            zb_in = Zbi
            sb_in = StBi

            idx = (SBT["Xb"] <= xb_in).sum()
            qb_in = SBT["Q"].iloc[idx - 1] if idx > 0 else np.nan

            zs_in = Zsi
            ss_in = StSi
            qs_in = sur_q(lb, alpha, ks, zs_in)

            # Storage Capacity Available
            ssa = sur_store(lb, alpha, ws, por, ws * alpha) - ss_in
            sba = max(SBT['S']) - sb_in  # Base Zone
            rech_in = min(recharge(lb, xb_in, ws, kz, zs_in, por), sba + qb_in)  # Initial Recharge Limited to Available Base Storage Capacity + Base Flow
            infil_in = 0

        # Initialize Time Step Using State Variables for Previous Time Step if Available
        if not ts_ini:
            xb_in = X[ts - 1]
            zb_in = Z[ts - 1, 1]
            sb_in = ST[ts - 1, 1]

            idx = (SBT["Xb"] <= xb_in).sum()
            qb_in = SBT["Q"].iloc[idx - 1] if idx > 0 else np.nan

            zs_in = Z[ts - 1, 0]
            ss_in = ST[ts - 1, 0]
            qs_in = sur_q(lb, alpha, ks, zs_in)

            # Storage Capacity Available
            ssa = sur_store(lb, alpha, ws, por, ws * alpha) - ss_in
            sba = max(SBT['S']) - sb_in  # Base Zone
            rech_in = min(recharge(lb, xb_in, ws, kz, zs_in, por), sba + qb_in)  # Initial Recharge Limited to Available Base Storage Capacity + Base Flow
            infil_in = 0


        I[ts] = 0
        # End of Time Step Calculations
        ss_en = max(ss_in + infil_in - rech_in - qs_in, 0)
        zs_en = sur_z(lb, alpha, ws, por, ss_en)
        qs_en = sur_q(lb, alpha, ks, zs_en)

        infil_en = min(infiltration(lb, ws, ks, alpha, zs_en, I[ts]), ssa)
        rech_en = min(recharge(lb, xb_in, ws, kz, zs_en, por), sba + qb_in)
        sb_en = max(sb_in + rech_en - qb_in, 0)
        idx = max((SBT["S"] < sb_en).sum(), 1) - 1

        # Safely extract the values from the DataFrame
        xb_en = SBT["Xb"].iloc[idx] if 0 <= idx < len(SBT) else np.nan
        zb_en = SBT["Z"].iloc[idx] if 0 <= idx < len(SBT) else np.nan
        qb_en = SBT["Q"].iloc[idx] if 0 <= idx < len(SBT) else np.nan

        # Final Calculations for Time Step
        if ts_ini:
            qcomp[ts, 0] = Surflow
            qcomp[ts, 1] = Baseflow
            EXC[ts, 1] = Rech

        if not ts_ini:
            qcomp[ts, 0] = (qs_in + qs_en) / 2  # Surface Flow
            qcomp[ts, 1] = (qb_in + qb_en) / 2  # Base Flow

            EXC[ts, 0] = (infil_in + infil_en) / 2
            EXC[ts, 1] = (rech_in + rech_en) / 2

        if ts_ini:
            ST[ts, 1] = StBi
            ST[ts, 0] = StSi
            Z[ts, 1] = Zbi
            Z[ts, 0] = Zsi

        # For Time Steps When States Are Available for Previous Time Step
        if not ts_ini:
            ST[ts, 0] = max(ST[ts - 1, 0] + EXC[ts, 0] - qcomp[ts, 0] - EXC[ts, 1], 0)
            ST[ts, 0] = min(ST[ts, 0], sur_store(lb, alpha, ws, por, ws * alpha))
            Z[ts, 0] = sur_z(lb, alpha, ws, por, ST[ts, 0])

            ST[ts, 1] = max(ST[ts - 1, 1] + EXC[ts, 1] - qcomp[ts, 1], 0)
            ST[ts, 1] = min(ST[ts, 1], max(SBT['S']))

            idx = max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1  # Adjust for 0-based indexing
            Z[ts, 1] = SBT['Z'].iloc[idx] if 0 <= idx < len(SBT) else np.nan

        #ETA[ts] = qin[ts] - np.sum(qcomp[ts, 0:3])  # Streamflow Residual
        idx = max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1  # Adjust for 0-based indexing
        X[ts] = SBT['Xb'].iloc[idx] if 0 <= idx < len(SBT) else np.nan
        ts += 1
        ts_ini = False
        #CLOSE CONDITION ts<p

    # OUTPUT
    tmp = pd.DataFrame({'Date': date, 'Baseflow': qcomp[:, 1],  'StSur': ST[:, 0], 'StBase': ST[:, 1], 'Zs': Z[:, 0], 'Zb': Z[:, 1], 'Rech': EXC[:, 1]})
    tmp = tmp[['Date', 'Baseflow', 'StSur', 'StBase', 'Zs', 'Zb', 'Rech']]
    return tmp


def plot_baseflow_simulation(streamflow, tmp, title="Baseflow Simulation"):
    """Plots observed streamflow vs simulated baseflow from PyBFS

    Creates a time series plot comparing observed total streamflow with simulated
    baseflow. Useful for visualizing baseflow separation results and assessing
    model performance.

    Parameters
    ----------
    streamflow : pd.DataFrame
        DataFrame containing observed data with columns Date (datetime) and Streamflow
        (observed streamflow m³/day)
    tmp : pd.DataFrame
        Output from PyBFS() containing flow components with column Baseflow (simulated
        baseflow component m³/day)
    title : str, optional
        Plot title, default is "Baseflow Simulation"

    Returns
    -------
    pd.DataFrame
        DataFrame used for plotting with columns date (datetime), streamflow (observed
        streamflow m³/s converted from m³/day), and baseflow (simulated baseflow m³/s
        converted from m³/day)

    Notes
    -----
    Streamflow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Black line shows total observed streamflow. Green line shows simulated baseflow component.
    Displays plot using matplotlib.

    Examples
    --------
    >>> results = PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
    >>> plot_df = plot_baseflow_simulation(streamflow_data, results)
    """
    # Prepare DataFrame for plotting
    df = pd.DataFrame({
        "date": pd.to_datetime(streamflow["Date"]),
        "streamflow": streamflow["Streamflow"],
        "baseflow": tmp["Baseflow"]
    })

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot observed streamflow (converted from m³/day to m³/s)
    ax.plot(df["date"], df["streamflow"] / 86400, color="black",
            label="Streamflow", linewidth=1)

    # Plot simulated baseflow (converted similarly)
    ax.plot(df["date"], df["baseflow"] / 86400, color="green",
            label="PyBFS", linewidth=1.5)

    # Labels and formatting
    ax.set_xlabel("Date", fontsize=16)
    ax.set_ylabel("Flow (cms)", fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.legend(loc="upper right", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=14)

    # Date formatting
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    fig.autofmt_xdate()

    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

    return df


def plot_forecast_baseflow(forecast_data):
    """Plot baseflow forecast time series

    Creates a time series plot of forecasted baseflow. Shows how baseflow
    is expected to evolve over the forecast period based on drainage from
    storage reservoirs.

    Parameters
    ----------
    forecast_data : pd.DataFrame
        Forecast results from forecast() function with columns Date (datetime of forecast)
        and Baseflow (forecasted baseflow m³/day)

    Notes
    -----
    Baseflow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Green line shows forecasted baseflow. Displays plot using matplotlib. No observed
    streamflow is shown (forecast period has no observations).

    Examples
    --------
    >>> forecast_result = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
    >>> plot_forecast_baseflow(forecast_result)
    """
    fig, axs = plt.subplots(figsize=(9, 5))
    date = pd.to_datetime(forecast_data["Date"])
    axs.plot(date, forecast_data['Baseflow']/86400, color='green', label='PyBFS Baseflow', linewidth=1.5)

    # Add legend
    axs.legend(loc='upper right', fontsize=13)

    # Set title and axis labels
    axs.set_title(f"Baseflow Forecast", fontsize=18)
    axs.set_xlabel('Date', fontsize=16)
    axs.set_ylabel('Flow (cms)', fontsize=16)

    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()

    # Tick label font sizes
    axs.tick_params(axis='x', labelsize=14)
    axs.tick_params(axis='y', labelsize=14)

    plt.show()


def plot_forecast_baseflow_streamflow(forecast_data, streamflow):
    """Plot forecast baseflow with observed streamflow for comparison

    Creates a time series plot comparing forecasted baseflow against observed
    streamflow. Useful for validating forecast performance when observations
    are available for the forecast period.

    Parameters
    ----------
    forecast_data : pd.DataFrame
        Forecast results from forecast() function with columns Date (datetime of forecast)
        and Baseflow (forecasted baseflow m³/day)
    streamflow : pd.DataFrame
        Observed streamflow data for the forecast period with columns Date (datetime of
        observations) and Streamflow (observed streamflow m³/day)

    Notes
    -----
    Flow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Blue line shows observed total streamflow. Green line shows forecasted baseflow.
    Displays plot using matplotlib. Used for forecast validation when observations become available.

    Examples
    --------
    >>> forecast_result = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
    >>> # Get observed data for same period
    >>> obs_data = streamflow_data[(streamflow_data['Date'] >= '2018-10-01') &
    ...                             (streamflow_data['Date'] <= '2018-11-30')]
    >>> plot_forecast_baseflow_streamflow(forecast_result, obs_data)
    """
    fig, axs = plt.subplots(figsize=(9, 5))
    date = pd.to_datetime(forecast_data["Date"])

    axs.plot(date, streamflow['Streamflow']/86400, color='blue', label='USGS Streamflow', linewidth=1.5)
    axs.plot(date, forecast_data['Baseflow']/86400, color='green', label='PyBFS Baseflow', linewidth=1.5)

    # Add legend
    axs.legend(loc='upper right', fontsize=13)

    # Set title and axis labels
    axs.set_title(f"Baseflow Forecast", fontsize=18)
    axs.set_xlabel('Date', fontsize=16)
    axs.set_ylabel('Flow (cms)', fontsize=16)

    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()

    # Tick label font sizes
    axs.tick_params(axis='x', labelsize=14)
    axs.tick_params(axis='y', labelsize=14)
    #end

    plt.show()

