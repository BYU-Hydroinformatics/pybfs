# -*- coding: utf-8 -*-
"""Utility functions for PyBFS

Contains helper functions for baseflow separation calculations including
surface and subsurface reservoir calculations, parameter extraction, and
baseflow table generation.
"""

import numpy as np
import pandas as pd
import math

from statsmodels.regression.quantile_regression import QuantReg

from numba import jit


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


# JIT-compiled versions of utility functions for performance
@jit(nopython=True, cache=True)
def sur_z_jit(lb, a, ws, por, ss):
    """JIT-compiled version of sur_z"""
    a1 = 1 / (2 * a)
    b1 = -2 * ws
    c1 = ss / (lb * por)
    
    discriminant = b1 ** 2 - 4 * a1 * c1
    
    if discriminant < 0:
        return ws * a
    else:
        return (-b1 - math.sqrt(discriminant)) / (2 * a1)


@jit(nopython=True, cache=True)
def sur_store_jit(lb, a, ws, por, zs):
    """JIT-compiled version of sur_store"""
    z = min(ws * a, zs)
    result = lb * (2 * ws * zs - zs**2 / a) * por
    return result


@jit(nopython=True, cache=True)
def sur_q_jit(lb, a, ks, z):
    """JIT-compiled version of sur_q"""
    return 2 * lb * z * a * ks


@jit(nopython=True, cache=True)
def dir_q_jit(lb, a, z, i):
    """JIT-compiled version of dir_q"""
    return 2 * lb * z / a * i


@jit(nopython=True, cache=True)
def infiltration_jit(lb, ws, ks, a, zs, i):
    """JIT-compiled version of infiltration"""
    return 2 * lb * (ws - zs / a) * min(i, ks)


@jit(nopython=True, cache=True)
def recharge_jit(lb, xb, ws, kz, zs, por):
    """JIT-compiled version of recharge"""
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

    # Create log-spaced discharge values (matching R: 1000 values, then prepend 0)
    qq = np.logspace(tmp_range[0], tmp_range[1], num=1000)
    qq = np.concatenate([[0], qq])  # Prepend 0 to match R's c(0,qq)

    # Calculate z based on the value of b
    if b != 0.5:
        z = ((qq * x1 * (2 * b - 1)) / (wb * kb * b ** 2)) ** (b / (2 * b - 1))
    else:
        z = np.exp(qq * 2 * x1 / (wb * kb))

    # Compute x and storage s
    x = x1 * z ** (1 / b)
    s = wb * por * ((1 / x1 ** b) * (1 / (b + 1)) * x ** (b + 1) + (lb - x) * z)

    # Construct the table and round values to 5 significant digits (matching R's signif(x,5))
    def signif(x, digits=5):
        """Round to specified number of significant digits (matching R's signif function)"""
        # Handle zero and NaN
        mask = (x == 0) | np.isnan(x)
        result = np.zeros_like(x)
        # For non-zero values, round to significant digits
        non_zero = ~mask
        if np.any(non_zero):
            # Calculate order of magnitude
            order = np.floor(np.log10(np.abs(x[non_zero])))
            # Round to significant digits
            result[non_zero] = np.round(x[non_zero] / (10 ** order), digits - 1) * (10 ** order)
        result[mask] = x[mask]  # Keep zeros and NaN as-is
        return result
    
    BT = pd.DataFrame({
        'Xb': signif(x, 5),
        'Z':  signif(z, 5),
        'S':  signif(s, 5),
        'Q':  signif(qq, 5)
    })

    # Filter out rows where x < lb
    BT = BT[BT['Xb'] < lb]

    return BT


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
    >>> result = bfs(streamflow_cal, SBT, basin_char, gw_hyd, flow)
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


def bf_ci(bfs_out):
    """Calculate credible intervals for baseflow separation

    Calculates credible intervals (5% to 95%) for streamflow given baseflow
    from BFS output. Used for uncertainty quantification in baseflow estimates.

    Parameters
    ----------
    bfs_out : pd.DataFrame
        Output DataFrame from bfs() function with columns:
        - Qsim.L3 or Qsim: Simulated total streamflow
        - Qob.L3 or Qob: Observed streamflow
        - DirectRunoff.L3 or DirectRunoff: Direct runoff component

    Returns
    -------
    tuple
        (ci_table, ci) where:
        - ci_table: DataFrame with credible intervals by flow quantile
        - ci: DataFrame with daily credible intervals

    Notes
    -----
    The function calculates fractional errors and groups them by flow quantiles
    to estimate uncertainty in baseflow separation.
    """
    # Handle column name variations
    qsim_col = 'Qsim' if 'Qsim' in bfs_out.columns else 'Qsim.L3'
    qob_col = 'Qob' if 'Qob' in bfs_out.columns else 'Qob.L3'
    dr_col = 'DirectRunoff' if 'DirectRunoff' in bfs_out.columns else 'DirectRunoff.L3'

    tmp_error = (bfs_out[qsim_col] - bfs_out[qob_col]) / bfs_out[qsim_col]
    tmp_error[np.abs(tmp_error) == np.inf] = 1
    tmp_error[bfs_out[dr_col] > 0] = np.nan

    tmp_q = bfs_out[qsim_col].copy()
    tmp_q[bfs_out[dr_col] > 0] = np.nan

    # Quantiles of simulated flow
    q_positive = tmp_q[tmp_q > 0]
    if len(q_positive) > 0:
        qnts = np.quantile(q_positive, np.arange(0.05, 1.0, 0.05))
    else:
        qnts = np.array([])

    # Standalone table
    ci_table = []
    for x in range(1, 18):  # x from 2 to 18 in R (0-indexed: 1 to 17)
        if x < len(qnts) - 1:
            y = (tmp_q > qnts[x-1]) & (tmp_q < qnts[x+1])
            error_subset = tmp_error[y]
            # Check if there are valid (non-NaN) values before computing quantiles
            if np.any(~np.isnan(error_subset)):
                ci_table.append([
                    qnts[x],
                    np.nanquantile(error_subset, 0.05),
                    np.nanquantile(error_subset, 0.50),
                    np.nanquantile(error_subset, 0.95)
                ])

    if ci_table:
        ci_table_df = pd.DataFrame(ci_table, columns=['Qsim.L3.T', 'FrEr0.05', 'FrEr0.50', 'FrEr0.95'])
    else:
        ci_table_df = pd.DataFrame(columns=['Qsim.L3.T', 'FrEr0.05', 'FrEr0.50', 'FrEr0.95'])

    # CI for each daily value
    qnt = np.full(len(tmp_q), np.nan)
    if len(qnts) >= 19 and len(ci_table_df) > 0:
        qnts_mid = (qnts[:18] + qnts[1:19]) / 2
        for t in range(len(tmp_q)):
            if bfs_out[dr_col].iloc[t] == 0:
                matches = bfs_out[qsim_col].iloc[t] < qnts_mid
                if np.any(matches):
                    qnt[t] = np.where(matches)[0][0]  # 0-indexed for Python

    # Create CI DataFrame
    ci = pd.DataFrame(index=bfs_out.index)
    if len(ci_table_df) > 0:
        ci['CB0.05'] = np.nan
        ci['CB0.95'] = np.nan
        for t in range(len(tmp_q)):
            if not np.isnan(qnt[t]) and int(qnt[t]) < len(ci_table_df):
                idx = int(qnt[t])
                ci.loc[ci.index[t], 'CB0.05'] = bfs_out[qsim_col].iloc[t] * (1 - ci_table_df['FrEr0.95'].iloc[idx])
                ci.loc[ci.index[t], 'CB0.95'] = bfs_out[qsim_col].iloc[t] * (1 - ci_table_df['FrEr0.05'].iloc[idx])
    else:
        ci['CB0.05'] = np.nan
        ci['CB0.95'] = np.nan

    return (ci_table_df, ci)


def flow_metrics(qin, timestep='day', fr4rise=0.05):
    """Calculate streamflow metrics used for baseflow separation

    Calculates flow metrics from a streamflow time series including recession
    coefficients, precision, and flow thresholds. These metrics are used to
    parameterize the baseflow separation algorithm.

    Parameters
    ----------
    qin : array-like
        Numeric vector with streamflow time series (m³/day or m³/hour)
    timestep : str, optional
        Either 'day' or 'hour', default is 'day'
    fr4rise : float, optional
        Fraction for rise detection threshold, default is 0.05

    Returns
    -------
    list
        Flow metrics [Qthresh, Rs, Rb1, Rb2, Prec, Fr4Rise] where:

        - Qthresh: Lower threshold on baseflow (m³/day or m³/hour)
        - Rs: Exponential coefficient for rapid recession (storm flow) [1/timestep]
        - Rb1: Exponential coefficient for slow recession (base flow at mean) [1/timestep]
        - Rb2: Exponential coefficient for intermediate recession (high base flow) [1/timestep]
        - Prec: Precision based on measurement resolution (m³/day or m³/hour)
        - Fr4Rise: Fraction for rise detection (same as input)

    Notes
    -----
    This function uses `statsmodels` for quantile regression.

    The function:
    1. Calculates precision from measurement resolution
    2. Interpolates missing values for gaps < 2 days
    3. Identifies recession periods
    4. Calculates 2-day recession rates (Rs)
    5. Calculates 10-day recession rates (Rb1, Rb2)
    6. Determines Qthresh from recession rate analysis

    Examples
    --------
    >>> import numpy as np
    >>> streamflow = np.array([10, 12, 15, 18, 16, 14, 12, 11, 10, 9])
    >>> metrics = flow_metrics(streamflow, timestep='day', fr4rise=0.05)
    >>> qthresh, rs, rb1, rb2, prec, fr4rise = metrics
    """
    qin = np.array(qin, dtype=float)
    xx = len(qin)

    # Set negative flows to NaN
    qin[qin < 0] = np.nan

    # PRECISION OF MEASURED BASEFLOW
    if np.any(~np.isnan(qin)):
        Q01 = np.nanquantile(qin, 0.01)
    else:
        Q01 = np.nan
    tmp = qin - Q01
    tmp_positive = tmp[tmp > 0]
    if len(tmp_positive) > 0:
        prec = np.nanmin(tmp_positive)  # Smallest positive difference
    else:
        prec = np.nan
    prec = max(prec, 0.1 * Q01)  # Increase if 10% of Q01 is larger
    
    # Debug output (commented out for production)
    # print(f"    Prec calculation debug:")
    # print(f"      Q01 (1st percentile): {Q01}")
    # print(f"      Smallest positive (qin - Q01): {np.nanmin(tmp_positive) if len(tmp_positive) > 0 else 'N/A'}")
    # print(f"      10% of Q01: {0.1 * Q01}")
    # print(f"      Final prec: {prec}")

    # Increase precision to median daily range for hourly data
    if timestep == 'hour':
        tmp_range = []
        for t in range(xx):
            start_idx = max(0, t - 23)
            end_idx = min(xx, t + 1)
            if end_idx > start_idx:
                daily_range = np.nanmax(qin[start_idx:end_idx]) - np.nanmin(qin[start_idx:end_idx])
                tmp_range.append(daily_range)
        if tmp_range:
            prec = max(prec, np.nanmedian(tmp_range))

    # Set window size for recession analysis
    if timestep == 'day':
        x = 2
    elif timestep == 'hour':
        x = 48
    else:
        raise ValueError("timestep must be 'day' or 'hour'")

    # INTERPOLATE STREAMFLOW IF MISSING PERIOD IS LESS THAN 2 DAYS
    tmp1 = np.arange(1, xx + 1)
    # Forward fill distance
    tmp2 = tmp1 - np.maximum.accumulate(tmp1 * (~np.isnan(qin)))
    # Backward fill distance
    tmp3_rev = tmp1[::-1] - np.maximum.accumulate(tmp1[::-1] * (~np.isnan(qin[::-1])))
    tmp3 = tmp3_rev[::-1]
    tmp = np.maximum(tmp2, tmp3)

    # Forward and backward fill
    qint_forward = qin.copy()
    qint_backward = qin.copy()
    for y in range(1, xx - 1):
        if np.isnan(qint_forward[y]):
            qint_forward[y] = qint_forward[y - 1]
        if np.isnan(qint_backward[xx - y - 1]):
            qint_backward[xx - y - 1] = qint_backward[xx - y]

    # Interpolate where gap < x
    mask = np.isnan(qin) & (tmp < x)
    qin[mask] = (qint_forward[mask] + qint_backward[mask]) / 2

    # IDENTIFY RECESSION DAYS - USING WINDOW OF PREVIOUS X TIME STEPS
    # R code: rp_indices=array(c(1:xx),dim=c(xx,1))
    #         for(y in 1:x){rp_indices=cbind(rp_indices,c((1-y):(xx-y)))}
    # R is 1-indexed: row y (from x+1 to xx) should have [y, y-1, ..., y-x]
    # For example, with x=2, row 3 should have [3, 2, 1]
    # Python 0-indexed: row y (from x to xx-1) should have [y, y-1, ..., y-x]
    # For example, with x=2, row 2 should have [2, 1, 0]
    rp_indices = np.arange(xx).reshape(-1, 1)  # First column: [0, 1, ..., xx-1] (current time step)
    for y_offset in range(1, x + 1):
        # R: c((1-y):(xx-y)) creates a vector where element at position i is (i+1-y)
        # For y=1: positions 0 to xx-2 get values 0 to xx-2
        # For y=2: positions 0 to xx-2 get values -1 to xx-3
        # But we need: for row i, column y_offset should be i - y_offset
        # So: column y_offset at row i should be i - y_offset
        # This means: indices[i] = i - y_offset for i from 0 to xx-1
        indices = np.arange(xx) - y_offset  # Each row i gets value i - y_offset
        rp_indices = np.column_stack([rp_indices, indices])

    rec = np.zeros(xx, dtype=bool)  # Logical vector indicating recession time steps
    rise = np.zeros(xx, dtype=bool)  # Logical vector indicating rise time steps
    qrise = np.zeros(xx)  # Rise over window
    qmax = np.zeros(xx, dtype=bool)  # Time step has max Q for window

    # Calculate change in streamflow over previous x time steps
    # R code: for(y in (x+1):xx) {qp=qin[rp_indices[y,]]; dq=qp[1:x]-qp[2:(x+1)]}
    # R is 1-indexed: rp_indices[y,] = [y, y-1, ..., y-x] for y from x+1 to xx
    # So qp = [qin[y], qin[y-1], ..., qin[y-x]]
    # dq = qp[1:x] - qp[2:(x+1)] = [qp[1], ..., qp[x]] - [qp[2], ..., qp[x+1]]
    #    = [qin[y-1], ..., qin[y-x]] - [qin[y-2], ..., qin[y-x-1]]
    # Python 0-indexed: for y from x to xx-1, rp_indices[y] = [y, y-1, ..., y-x]
    # So qp = [qin[y], qin[y-1], ..., qin[y-x]]
    # dq should be: [qp[0], ..., qp[x-1]] - [qp[1], ..., qp[x]]
    #    = [qin[y], ..., qin[y-x+1]] - [qin[y-1], ..., qin[y-x]]
    # Wait, that's not right. Let me match R exactly:
    # R: dq = qp[1:x] - qp[2:(x+1)] where qp indices are 1-indexed
    # Python: dq = qp[1:x+1] - qp[2:x+2] where qp indices are 0-indexed
    # Actually, in R, qp[1:x] means elements 1 through x (1-indexed)
    # In Python, that's qp[0:x] (0-indexed)
    # So: dq = qp[0:x] - qp[1:x+1]
    for y in range(x, xx):
        idx_row = rp_indices[y, :x+1]  # Get first x+1 columns (current + x previous)
        valid_idx = idx_row[idx_row >= 0]  # Filter valid indices
        if len(valid_idx) == x + 1:
            qp = qin[valid_idx]
            if not np.any(np.isnan(qp)):
                # R: dq = qp[1:x] - qp[2:(x+1)] (1-indexed)
                # Python: dq = qp[0:x] - qp[1:x+1] (0-indexed)
                # This gives: dq = [qp[0]-qp[1], qp[1]-qp[2], ..., qp[x-1]-qp[x]]
                # Which is: [qin[y]-qin[y-1], qin[y-1]-qin[y-2], ..., qin[y-x+1]-qin[y-x]]
                dq = qp[0:x] - qp[1:x+1]  # Change for each time step in window
                q4rise = max(fr4rise * qp[0], prec)  # Threshold for rise
                qrise[y] = np.nansum(dq)  # Increase for window
                rise[y] = (qrise[y] > q4rise)
                if qrise[y] < q4rise:
                    qrise[y] = np.nan

                # Recession: all changes < prec and cumulative change < 0
                rec_cond1 = np.all(dq < prec)
                rec_cond2 = (np.nansum(dq) < 0)
                rec[y] = rec_cond1 & rec_cond2
                qmax[y] = (qp[0] == np.nanmax(qp))  # Time step has maximum streamflow
                
                # Debug output (commented out for production)
                # if y == x:  # First valid period
                #     print(f"      Debug first period y={y}:")
                #     print(f"        idx_row={idx_row}")
                #     print(f"        valid_idx={valid_idx}")
                #     print(f"        qp={qp}")
                #     print(f"        dq={dq}")
                #     print(f"        qp[0:x]={qp[0:x]}, qp[1:x+1]={qp[1:x+1]}")

    # 2-DAY RECESSION RATES FOLLOWING A RISE
    # R code: dq=c(rep(NA,x),(qin[(x+1):xx]-qin[1:(xx-x)])/qrise[1:(xx-x)])
    # R is 1-indexed: for k from 1 to (xx-x):
    #   dq[x+k] = (qin[x+k] - qin[k]) / qrise[k]
    # Note: qrise[k] for k=1..x is 0 (unset), for k=(x+1)..(xx-x) is set
    # Python 0-indexed: for k from 0 to (xx-x-1):
    #   dq[x+k] = (qin[x+k] - qin[k]) / qrise[k]
    # Note: qrise[k] for k=0..(x-1) is 0 (unset), for k=x..(xx-x-1) is set
    dq = np.full(xx, np.nan)
    # Calculate for indices from x to xx-1
    for k in range(xx - x):
        i = x + k  # Current index in dq and qin (Python 0-indexed)
        # Match R exactly: use qrise[k] even if it's 0 (will produce Inf, filtered later)
        if k < len(qrise):
            with np.errstate(divide='ignore', invalid='ignore'):
                dq[i] = (qin[i] - qin[k]) / qrise[k]
    
    dq[~rec] = np.nan  # Limit to recession periods
    # Set to NaN for time step with max streamflow at start of window
    # R code: dq[!c(rep(FALSE,x),qmax[1:(xx-x)])]=NA
    # qmax is only set for indices (x+1):xx, so qmax[1:(xx-x)] includes:
    #   - qmax[1] through qmax[x]: FALSE (unset)
    #   - qmax[x+1] through qmax[xx-x]: set values
    # In Python: qmax is set for indices x:xx-1, so we need qmax[0:(xx-x)]
    # But qmax[0] through qmax[x-1] are False (unset), matching R's qmax[1] through qmax[x]
    qmax_mask = np.concatenate([np.zeros(x, dtype=bool), qmax[:xx-x]])
    dq[~qmax_mask] = np.nan

    # Check if there are valid values before computing quantile
    dq_valid = dq[~np.isnan(dq)]
    num_recession = np.sum(rec)
    num_qmax = np.sum(qmax)
    num_valid_dq = len(dq_valid)
    
    # Debug output (commented out for production)
    # print(f"    Rs calculation debug:")
    # print(f"      Total time steps: {xx}")
    # print(f"      Recession periods identified: {num_recession}")
    # print(f"      Max Q periods: {num_qmax}")
    # print(f"      Valid dq values after filtering: {num_valid_dq}")
    # print(f"      prec value used: {prec:.2f}")
    
    if len(dq_valid) > 0:
        dq_95 = np.nanquantile(dq, 0.95)
    else:
        dq_95 = np.nan
    
    if np.isnan(dq_95):
        # If no valid recession rates, use default
        if num_recession > 0 and num_valid_dq == 0:
            Rs = -0.05  # Less aggressive default
        else:
            Rs = -0.1  # Default when no recession data
    else:
        if dq_95 < -1:
            dq_95 = -0.9  # Reduce 100% recession to 90%
        if dq_95 <= -1:
            dq_95 = -0.9  # Ensure we can take log(1 + dq_95)
        Rs = np.log(1 + dq_95) / x  # Recession coefficient for exponential model
        if np.isnan(Rs) or np.isinf(Rs):
            Rs = -0.1  # Default value

    # 10-DAY RECESSIONS
    if timestep == 'day':
        x = 10
    elif timestep == 'hour':
        x = 240

    rec = np.zeros(xx, dtype=bool)  # Reset recession vector

    # Rebuild indices for 10-day window (0-based for Python)
    # R: rp_indices=array(c(1:xx),dim=c(xx,1)); for(y in 1:x){rp_indices=cbind(rp_indices,c((1-y):(xx-y)))}
    # R is 1-indexed: first column is [1, 2, ..., xx]
    # Python 0-indexed: first column should be [0, 1, ..., xx-1]
    rp_indices = np.arange(xx).reshape(-1, 1)  # [0, 1, ..., xx-1]
    for y in range(1, x + 1):
        # R: c((1-y):(xx-y))
        # For y=1: [0, 1, ..., xx-2] (1-indexed: 1-1=0 to xx-1)
        # Python: indices should be [0-y, 1-y, ..., xx-1-y] = [-y, -y+1, ..., xx-1-y]
        indices = np.arange(xx) - y
        rp_indices = np.column_stack([rp_indices, indices])

    # Calculate change over previous x time steps
    # R: for(y in (x+1):xx){qp=qin[rp_indices[y,]]; if(any(!is.na(qp))){dq=qp[1:x]-qp[2:(x+1)]; ...; rec[y]=(qp[1]==min(qp,na.rm=T))}}
    # R is 1-indexed: for y from x+1 to xx
    # Python 0-indexed: for y from x to xx-1
    # rp_indices[y] in Python (0-indexed) = rp_indices[y+1,] in R (1-indexed)
    for y in range(x, xx):
        idx_row = rp_indices[y]  # Get row y (0-indexed), which is row y+1 in R (1-indexed)
        # Filter out negative indices (invalid)
        valid_mask = idx_row >= 0
        if np.any(valid_mask):
            qp = qin[idx_row[valid_mask]]
            if np.any(~np.isnan(qp)):
                # R: rec[y]=(qp[1]==min(qp,na.rm=T))
                # In R, qp[1] is the second element (1-indexed), which is qp[0] in Python (0-indexed)
                # But wait, rp_indices[y,] in R gives [y, y-1, ..., y-x] (1-indexed)
                # So qp = [qin[y], qin[y-1], ..., qin[y-x]] (1-indexed)
                # qp[1] = qin[y-1] (1-indexed)
                # Actually, let me check: rp_indices is built with first column [1,2,...,xx]
                # and subsequent columns are offsets. So rp_indices[y,] = [y, y-1, ..., y-x]
                # But the first column is just the row number, so rp_indices[y,1] = y
                # So qp[1] in R refers to the second column, which is y-1
                # But we want the "current" time step, which should be y
                # I think qp[1] in R actually refers to the element at index 1 in the qp array
                # which is qin[rp_indices[y,2]] = qin[y-1]
                # But the comment says "qp[1]" is the current time step...
                # Let me trust the R code: rec[y]=(qp[1]==min(qp,na.rm=T))
                # This means we check if the second element of qp equals the minimum
                # Since qp = qin[rp_indices[y,]] and rp_indices[y,] = [y, y-1, ..., y-x]
                # qp[1] = qin[y-1] (1-indexed) = qin[y-1] (0-indexed)
                # But that doesn't make sense as the "current" time step
                # Actually, I think the issue is that rp_indices includes the current time step as the first element
                # So qp[1] is actually the previous time step, not the current
                # But the comment says "TIME STEP MUST HAVE MINIMUM STREAMFLOW FOR WINDOW"
                # So maybe qp[1] is being compared to min(qp) to check if it's the minimum
                # Let me just match R exactly: use qp[0] (first element) as the current time step
                # and check if it equals the minimum
                # Actually, wait - in R, arrays are 1-indexed, so qp[1] is the first element
                # So qp[1] = qin[rp_indices[y,1]] = qin[y] (1-indexed) = qin[y] (0-indexed)
                # So in Python, we should use qp[0] and check if it equals min(qp)
                if len(qp) > 0:
                    rec[y] = (qp[0] == np.nanmin(qp))

    # 10-DAY RECESSION RATES
    # R: dq=c(rep(NA,x),(qin[(x+1):xx]-qin[1:(xx-x)])/qin[1:(xx-x)])
    # R is 1-indexed: dq[x+1:xx] = (qin[x+1:xx] - qin[1:(xx-x)]) / qin[1:(xx-x)]
    # Python 0-indexed: dq[x:xx] = (qin[x:xx] - qin[0:(xx-x)]) / qin[0:(xx-x)]
    dq = np.full(xx, np.nan)
    # Calculate for indices from x to xx-1
    for i in range(x, xx):
        if i - x < len(qin) and qin[i - x] > 0 and not np.isnan(qin[i - x]):
            dq[i] = (qin[i] - qin[i - x]) / qin[i - x]
    dq[~rec] = np.nan  # Limit to recession periods

    # Limit to positive streamflow where stream did not dry completely (dq > -1)
    # R: tmp.d=qin>0; tmp.d[is.na(qin)]=FALSE; tmp.d[dq==-1]=FALSE
    # Note: R excludes dq==-1 (exact equality), not dq<=-1
    tmp_d = (qin > 0) & (~np.isnan(qin)) & (~np.isnan(dq)) & (dq != -1)
    if np.any(tmp_d):
        logq = np.log10(qin[tmp_d])
        rb = np.log10(1 + dq[tmp_d]) / x
    else:
        # No valid data for quantile regression
        logq = np.array([])
        rb = np.array([])

    # SET Qthresh TO MINIMUM STREAMFLOW EXCLUDING NO FLOW PERIODS
    q_positive = qin[qin > 0]
    if len(q_positive) == 0:
        # No positive flows - return default values
        return [np.nan, np.nan, np.nan, np.nan, np.nan, fr4rise, np.array([np.nan, np.nan])]
    
    Qthresh = np.nanmin(q_positive)
    qmean = np.nanmean(qin)
    
    if np.isnan(Qthresh) or np.isnan(qmean) or Qthresh <= 0 or qmean <= 0:
        return [np.nan, np.nan, np.nan, np.nan, np.nan, fr4rise, np.array([np.nan, np.nan])]

    # QUANTILES OF STREAMFLOW, EXCLUDING NO-FLOW PERIODS, UP TO MEAN
    q_positive = qin[(qin > 0) & (qin <= qmean)]
    if len(q_positive) > 0:
        tmp_qq = np.quantile(q_positive, np.arange(0, 1.01, 0.01))

        # Calculate median recession rates for streamflow quantiles
        # R: for(qq in tmp.qq) {tmp.r=c(tmp.r,quantile(dq[qin<qq],p=0.95,na.rm=T))}
        # Note: dq should only contain recession rates (negative values), but we filter to be safe
        tmp_r = []
        for qq in tmp_qq:
            # Filter: qin < qq and dq is valid and negative (recession)
            recession_at_qq = dq[(qin < qq) & (~np.isnan(dq)) & (~np.isinf(dq)) & (dq < 0)]
            if len(recession_at_qq) > 0 and np.any(~np.isnan(recession_at_qq)):
                tmp_r.append(np.nanquantile(recession_at_qq, 0.95))
            else:
                tmp_r.append(np.nan)

        # Set Qthresh to streamflow with maximum (least negative) recession rate
        if len(tmp_r) > 0 and not np.all(np.isnan(tmp_r)):
            tmp_qqthresh_idx = np.nanargmax(tmp_r)
            tmp_qthresh = tmp_qq[tmp_qqthresh_idx]
            Qthresh = max(Qthresh, tmp_qthresh)

    # LOWER LIMIT ON FASTEST RECESSION RATES AS A FUNCTION OF LOG STREAMFLOW
    # Filter out NaN and Inf values
    if len(logq) == 0 or len(rb) == 0:
        rb10 = np.array([np.nan, np.nan])
    else:
        valid_mask = ~(np.isnan(logq) | np.isnan(rb) | np.isinf(logq) | np.isinf(rb))
        logq_clean = logq[valid_mask]
        rb_clean = rb[valid_mask]

        if len(rb_clean) < 9 or len(logq_clean) < 9:
            rb10 = np.array([np.nan, np.nan])
        elif np.nanmin(logq_clean) == np.nanmax(logq_clean):
            rb10 = np.array([np.nan, np.nan])
        else:
            try:
                # Quantile regression at 10th percentile (tau=0.1)
                # rb ~ logq: rb is dependent (y), logq is independent (X)
                # Add intercept term explicitly (column of ones)
                X = np.column_stack([np.ones(len(logq_clean)), logq_clean.reshape(-1, 1)])
                y = rb_clean
                model = QuantReg(y, X)
                result = model.fit(q=0.1)
                # Check if we have both intercept and slope
                if len(result.params) >= 2:
                    rb10 = np.array([result.params[0], result.params[1]])  # intercept, slope
                elif len(result.params) == 1:
                    # Only intercept, no slope - use default slope
                    rb10 = np.array([result.params[0], 0.0])
                else:
                    rb10 = np.array([np.nan, np.nan])
            except Exception:
                # If quantile regression fails, return NaN
                rb10 = np.array([np.nan, np.nan])

    # RECESSION RATE FOR MEAN STREAMFLOW
    if not np.any(np.isnan(rb10)) and not np.any(np.isinf(rb10)):
        try:
            Rb1 = rb10[0] + rb10[1] * np.log10(qmean)
            if np.isnan(Rb1) or np.isinf(Rb1):
                Rb1 = -0.002
            else:
                Rb1 = min(Rb1, -0.002)  # Limit to -0.002
        except:
            Rb1 = -0.002  # Default if calculation failed
    else:
        Rb1 = -0.002  # Default if regression failed

    # RECESSION RATE FOR LOW BASEFLOW (10TH PERCENTILE OF FILTERED STREAMFLOW)
    if not np.any(np.isnan(rb10)) and not np.any(np.isinf(rb10)) and len(logq_clean) > 0:
        try:
            # logq_clean is already filtered, but use nanquantile to match R behavior
            tmp_logq10 = np.nanquantile(logq_clean, 0.1)
            if not np.isnan(tmp_logq10):
                Rb2 = rb10[0] + rb10[1] * tmp_logq10
                if np.isnan(Rb2) or np.isinf(Rb2):
                    Rb2 = -0.001
                else:
                    Rb2 = min(Rb2, -0.001)  # Limit to -0.001
            else:
                Rb2 = -0.001
        except:
            Rb2 = -0.001  # Default if calculation failed
    else:
        Rb2 = -0.001  # Default if regression failed

    # Validate all output values
    if np.any(np.isnan([Qthresh, Rs, Rb1, Rb2, prec])) or np.any(np.isinf([Qthresh, Rs, Rb1, Rb2, prec])):
        # Return NaN values if any are invalid
        return [np.nan, np.nan, np.nan, np.nan, np.nan, fr4rise, np.array([np.nan, np.nan])]
    
    # Return flow metrics and rb10
    # Note: rb10 is returned as 7th element for calibration functions
    return [Qthresh, Rs, Rb1, Rb2, prec, fr4rise, rb10]

