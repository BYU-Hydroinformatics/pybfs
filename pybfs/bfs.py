# -*- coding: utf-8 -*-
"""Baseflow Separation (BFS) function

Main function for performing physically-based baseflow separation using
a coupled surface-subsurface reservoir model.
"""

import numpy as np
import pandas as pd
import math
from .utilities import (
    sur_z, sur_store, sur_q, dir_q, infiltration, recharge
)

from numba import jit
from .utilities import (
    sur_z_jit, sur_store_jit, sur_q_jit, dir_q_jit,
    infiltration_jit, recharge_jit
)
HAS_NUMBA = True


if HAS_NUMBA:
    @jit(nopython=True, cache=True)
    def _bfs_core_loop(qin, sbt_xb, sbt_z, sbt_s, sbt_q, rise, p, 
                       lb, alpha, ws, por, ks, kz, prec, ifact, qmean,
                       X, qcomp, ETA, I, Z, ST, EXC):
        """
        JIT-compiled core computation loop for BFS algorithm.
        This is the computationally intensive part that benefits from JIT.
        """
        # Find first valid time step
        ts = 0
        while ts < p and math.isnan(qin[ts]):
            ts += 1
        
        if ts >= p:
            return
        
        # Initialization
        ts_ini = True
        qb_in = min(qin[ts], qmean)
        
        # Find index in SBT - match original: (SBT["Q"] <= qb_in).sum() - 1
        idx = 0
        for i in range(len(sbt_q)):
            if sbt_q[i] <= qb_in:
                idx = i + 1  # Count (like .sum())
            else:
                break
        idx = idx - 1 if idx > 0 else -1  # Convert count to index (idx - 1)
        
        xb_in = sbt_xb[idx] if idx >= 0 and idx < len(sbt_xb) else np.nan
        zb_in = sbt_z[idx] if idx >= 0 and idx < len(sbt_z) else np.nan
        sb_in = sbt_s[idx] if idx >= 0 and idx < len(sbt_s) else np.nan
        
        # Calculate available base storage
        sba = 0.0
        for i in range(len(sbt_s)):
            if sbt_s[i] > sba:
                sba = sbt_s[i]
        sba = sba - sb_in
        
        # Surface flow calculations
        qs_in_val = qin[ts] - qb_in
        if math.isnan(qs_in_val):
            qs_in = 0.0
        else:
            qs_in = max(0.0, qs_in_val)
        
        zs_in_val = qs_in / (2 * lb * ks * alpha)
        if math.isnan(zs_in_val):
            zs_in = ws * alpha
        else:
            zs_in = min(zs_in_val, ws * alpha)
        
        ss_in = sur_store_jit(lb, alpha, ws, por, zs_in)
        ssa = sur_store_jit(lb, alpha, ws, por, ws * alpha) - ss_in
        
        infil_in = 0.0
        rech_in = recharge_jit(lb, xb_in, ws, kz, zs_in, por)
        
        # Main time step loop - this is the bottleneck
        while ts < p:
            if not ts_ini:
                xb_in = X[ts - 1]
                zb_in = Z[ts - 1, 1]
                sb_in = ST[ts - 1, 1]
                
                # Find index in SBT - match original: (SBT["Xb"] <= xb_in).sum() - 1
                idx = 0
                for i in range(len(sbt_xb)):
                    if sbt_xb[i] <= xb_in:
                        idx = i + 1  # Count (like .sum())
                    else:
                        break
                idx = idx - 1 if idx > 0 else -1  # Convert count to index (idx - 1)
                qb_in = sbt_q[idx] if idx >= 0 and idx < len(sbt_q) else np.nan
                
                zs_in = Z[ts - 1, 0]
                ss_in = ST[ts - 1, 0]
                qs_in = sur_q_jit(lb, alpha, ks, zs_in)
                
                ssa = sur_store_jit(lb, alpha, ws, por, ws * alpha) - ss_in
                sba = 0.0
                for i in range(len(sbt_s)):
                    if sbt_s[i] > sba:
                        sba = sbt_s[i]
                sba = sba - sb_in
                rech_in = min(recharge_jit(lb, xb_in, ws, kz, zs_in, por), sba + qb_in)
                qd = 0.0
                infil_in = 0.0
                
                I[ts] = 0.0
                etaest = max(0.0, qin[ts] - qb_in - qs_in)
                
                # Impulse calculation loop - nested loops, computationally intensive
                if ts > 0 and etaest > 0:
                    if rise[ts] or (ts > 0 and rise[ts - 1]):
                        I[ts] = etaest / (2 * lb * ws)
                        zs = zs_in
                        qs = qs_in
                        
                        for x in ifact:
                            i = I[ts]
                            eta = etaest
                            while eta > max(prec, qin[ts] / 100) and i > 0:
                                I[ts] = i
                                etaest = eta
                                i = x * i
                                infil = min(infiltration_jit(lb, ws, ks, alpha, (zs_in + zs) / 2, i), ssa)
                                ss = max(ss_in + infil - rech_in - qs, 0.0)
                                zs = sur_z_jit(lb, alpha, ws, por, ss)
                                qs = sur_q_jit(lb, alpha, ks, zs)
                                qd = (dir_q_jit(lb, alpha, zs_in, i) + 
                                      dir_q_jit(lb, alpha, (zs - zs_in), i / 2) + 
                                      max(2 * lb * (ws - zs_in / alpha) * (I[ts] - ks), 0.0))
                                eta = qin[ts] - qs - qd - qb_in
                
                infil_in = min(infiltration_jit(lb, ws, ks, alpha, zs_in, I[ts]), ssa)
            
            # End of time step calculations
            ss_en = max(ss_in + infil_in - rech_in - qs_in, 0.0)
            zs_en = sur_z_jit(lb, alpha, ws, por, ss_en)
            qs_en = sur_q_jit(lb, alpha, ks, zs_en)
            
            infil_val = infiltration_jit(lb, ws, ks, alpha, zs_en, I[ts])
            if math.isnan(infil_val):
                infil_en = ssa
            else:
                infil_en = min(infil_val, ssa)
            
            rech_en = min(recharge_jit(lb, xb_in, ws, kz, zs_en, por), sba + qb_in)
            sb_en = max(sb_in + rech_en - qb_in, 0.0)
            
            # Find index in SBT - match original: max((SBT["S"] < sb_en).sum(), 1) - 1
            idx = 0
            for i in range(len(sbt_s)):
                if sbt_s[i] < sb_en:
                    idx = i + 1  # Count (like .sum())
                else:
                    break
            idx = max(idx, 1) - 1  # max(count, 1) - 1
            
            xb_en = sbt_xb[idx] if 0 <= idx < len(sbt_xb) else np.nan
            zb_en = sbt_z[idx] if 0 <= idx < len(sbt_z) else np.nan
            qb_en = sbt_q[idx] if 0 <= idx < len(sbt_q) else np.nan
            
            # Final calculations
            qcomp[ts, 0] = (qs_in + qs_en) / 2
            qcomp[ts, 1] = (qb_in + qb_en) / 2
            
            EXC[ts, 0] = (infil_in + infil_en) / 2
            EXC[ts, 1] = (rech_in + rech_en) / 2
            
            if ts_ini:
                ST[ts, 0] = ss_en
                ST[ts, 1] = sb_en
                Z[ts, 0] = zs_en
                Z[ts, 1] = zb_en
            
            if not ts_ini:
                ST[ts, 0] = max(ST[ts - 1, 0] + EXC[ts, 0] - qcomp[ts, 0] - EXC[ts, 1], 0.0)
                ST[ts, 0] = min(ST[ts, 0], sur_store_jit(lb, alpha, ws, por, ws * alpha))
                Z[ts, 0] = sur_z_jit(lb, alpha, ws, por, ST[ts, 0])
                ST[ts, 1] = max(ST[ts - 1, 1] + EXC[ts, 1] - qcomp[ts, 1], 0.0)
                
                max_s = 0.0
                for i in range(len(sbt_s)):
                    if sbt_s[i] > max_s:
                        max_s = sbt_s[i]
                ST[ts, 1] = min(ST[ts, 1], max_s)
                
                # Find index in SBT - match original: max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1
                idx = 0
                for i in range(len(sbt_s)):
                    if sbt_s[i] <= ST[ts, 1]:
                        idx = i + 1  # Count (like .sum())
                    else:
                        break
                idx = max(idx, 1) - 1  # max(count, 1) - 1
                Z[ts, 1] = sbt_z[idx] if idx >= 0 and idx < len(sbt_z) else np.nan
                
                qcomp[ts, 2] = (dir_q_jit(lb, alpha, zs_in, I[ts]) + 
                               dir_q_jit(lb, alpha, (Z[ts, 0] - zs_in), I[ts] / 2) + 
                               max(2 * lb * (ws - zs_in / alpha) * (I[ts] - ks), 0.0))
            
            ETA[ts] = qin[ts] - (qcomp[ts, 0] + qcomp[ts, 1] + qcomp[ts, 2])
            
            # Find index in SBT - match original: max((SBT['S'] <= ST[ts, 1]).sum(), 1) - 1
            idx = 0
            for i in range(len(sbt_s)):
                if sbt_s[i] <= ST[ts, 1]:
                    idx = i + 1  # Count (like .sum())
                else:
                    break
            idx = max(idx, 1) - 1  # max(count, 1) - 1
            X[ts] = sbt_xb[idx] if idx >= 0 and idx < len(sbt_xb) else np.nan
            
            ts += 1
            ts_ini = False


def bfs(streamflow, SBT, basin_char, gw_hyd, flow, timestep='day', error_basis='total', use_jit=True):
    """Main BFS function for baseflow separation

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
    timestep : str, optional
        Time step for calculations ('day' or 'hour'), default 'day'
    error_basis : str, optional
        Basis for error calculation ('base' or 'total'), default 'total'
    use_jit : bool, optional
        Whether to use NUMBA JIT compilation for faster computation (default True).
        If False, uses the non-JIT Python implementation.

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
    >>> results = bfs(streamflow_data, SBT, basin_char, gw_hyd, flow)
    """
    date = pd.to_datetime(streamflow["Date"])
    qin = np.array(streamflow['Streamflow'])
    # Remove negative flow values (matching R: qin[qin<0]=NA)
    qin[qin < 0] = np.nan
    qmean = np.nanmean(qin)

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
    # Handle divide by zero: when qin is 0, set dqfr based on conditions
    dqfr = np.zeros_like(dq, dtype=float)
    mask_zero_both = (dq == 0) & (qin == 0)
    mask_neg_dq_zero_qin = (dq < 0) & (qin == 0)
    mask_valid = qin != 0
    
    dqfr[mask_zero_both] = 0
    dqfr[mask_neg_dq_zero_qin] = 1
    dqfr[mask_valid] = dq[mask_valid] / qin[mask_valid]

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
    # Handle NaN values (matching R's na.rm=T)
    qs_in_val = qin[ts] - qb_in
    if np.isnan(qs_in_val):
        qs_in = 0
    else:
        qs_in = max(0, qs_in_val)  # Ensure surface flow is non-negative
    
    zs_in_val = qs_in / (2 * lb * ks * alpha)
    if np.isnan(zs_in_val):
        zs_in = ws * alpha
    else:
        zs_in = min(zs_in_val, ws * alpha)  # Saturated thickness of surface reservoir
    ss_in = sur_store(lb, alpha, ws, por, zs_in)  # Surface storage
    ssa = sur_store(lb, alpha, ws, por, ws * alpha) - ss_in  # AVAILABLE SURFACE STORAGE

    # Infiltration and recharge
    infil_in = 0
    rech_in = recharge(lb, xb_in, ws, kz, zs_in, por)

    # Use JIT-compiled core loop if requested
    if use_jit:
        # Convert SBT DataFrame to numpy arrays for JIT
        sbt_xb = np.array(SBT['Xb'].values)
        sbt_z = np.array(SBT['Z'].values)
        sbt_s = np.array(SBT['S'].values)
        sbt_q = np.array(SBT['Q'].values)
        ifact_arr = np.array(ifact)
        
        # Call JIT-compiled core loop
        _bfs_core_loop(qin, sbt_xb, sbt_z, sbt_s, sbt_q, rise, p,
                       lb, alpha, ws, por, ks, kz, prec, ifact_arr, qmean,
                       X, qcomp, ETA, I, Z, ST, EXC)
    else:
        # Original implementation - keep the while loop as-is
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
            # Handle NaN values in infil_en calculation (matching R's na.rm=T)
            infil_val = infiltration(lb, ws, ks, alpha, zs_en, I[ts])
            if np.isnan(infil_val):
                infil_en = ssa
            else:
                infil_en = min(infil_val, ssa)
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
                # Direct runoff is NOT set for initial time step in R (remains NA), so leave as NaN
                # qcomp[ts, 2] remains np.nan

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

    # Calculate Qsim - if any component is NaN, result should be NaN (matching R's rowSums behavior)
    qsim = np.sum(qcomp, axis=1)
    # R's rowSums returns NA if any element is NA, so set to NaN if any component is NaN
    qsim[np.isnan(qcomp).any(axis=1)] = np.nan
    
    tmp = pd.DataFrame({'Date': date, 'Qob': qin, 'Qsim': qsim, 'SurfaceFlow': qcomp[:, 0], 'Baseflow': qcomp[:, 1], 'DirectRunoff': qcomp[:, 2], 'X': X,'Eta': ETA, 'StSur': ST[:, 0], 'StBase': ST[:, 1], 'Impulse.L': I, 'Zs.L': Z[:, 0], 'Zb.L': Z[:, 1], 'Infil': EXC[:, 0], 'Rech': EXC[:, 1], 'RecessCount.T': recess_day, 'AdjPctEr': APE, 'Weight': Weight})
    tmp = tmp[['Date', 'Qob', 'Qsim', 'SurfaceFlow', 'Baseflow', 'DirectRunoff',  'X','Eta', 'StSur', 'StBase', 'Impulse.L', 'Zs.L', 'Zb.L', 'Infil', 'Rech', 'RecessCount.T', 'AdjPctEr', 'Weight']]
    return tmp

