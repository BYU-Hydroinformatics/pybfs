# -*- coding: utf-8 -*-
"""Skill assessment functions for PyBFS

Functions to evaluate how well pybfs performs baseflow separation and forecasting.
"""

import numpy as np
import pandas as pd
from .bfs import bfs
from .utilities import forecast


def modified_strict_baseflow(hydrograph, quantile=0.8):
    """Identify baseflow-dominated days using the modified strict baseflow method.

    A statistical method that uses streamflow derivatives and a quantile threshold
    to flag days where observed streamflow is likely composed entirely of baseflow.
    Does not use any pybfs model parameters.

    Parameters
    ----------
    hydrograph : pd.DataFrame
        Observed streamflow with columns 'Date' (datetime) and 'Streamflow' (m³/day).
    quantile : float, optional
        Quantile value for the strict baseflow filter (default 0.8).

    Returns
    -------
    np.ndarray
        Boolean array of length len(hydrograph). True = baseflow day, False = non-baseflow day.
    """
    Q = np.array(hydrograph['Streamflow'], dtype=float)
    n = len(Q)

    # Replace NaN with 0 temporarily for derivative calculations; mark NaN days as non-baseflow
    nan_mask = np.isnan(Q)
    Q_filled = Q.copy()
    Q_filled[nan_mask] = 0.0

    def _strict_baseflow(Q, quantile):
        dQ = (Q[2:] - Q[:-2]) / 2

        # 1. flow data associated with positive and zero values of dy/dt
        wet1 = np.concatenate([[True], dQ >= 0, [True]])

        # 2. previous 2 points before points with dy/dt >= 0, and next 3 points
        idx_first = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == 1)[0] + 1
        idx_last = np.where(wet1[1:].astype(int) - wet1[:-1].astype(int) == -1)[0]
        idx_before = np.repeat([idx_first], 2) - np.tile(range(1, 3), idx_first.shape)
        idx_next = np.repeat([idx_last], 3) + np.tile(range(1, 4), idx_last.shape)
        idx_remove = np.concatenate([idx_before, idx_next])
        wet2 = np.full(Q.shape, False)
        wet2[idx_remove.clip(min=0, max=Q.shape[0] - 1)] = True

        # 3. recession limbs in the upper quantile
        wet3_core = np.concatenate([[True], dQ[1:] - dQ[:-1] < 0, [True, True]])
        idx3_all = np.where(wet3_core)[0]
        q_thresh = np.quantile(Q, quantile)
        idx3 = idx3_all[Q[idx3_all] >= q_thresh]
        idx3_buffer = np.repeat([idx3], 10) + np.tile(range(-4, 6), idx3.shape)
        wet3 = np.full(Q.shape, False)
        wet3[idx3_buffer.clip(min=0, max=Q.shape[0] - 1)] = True

        # 4. flow data followed by a data point with a larger value of -dy/dt
        wet4 = np.concatenate([[True], dQ[1:] - dQ[:-1] < 0, [True, True]])

        dry = ~(wet1 | wet2 | wet3 | wet4)
        return dry

    def _flow_thresh(Q, strict_mask):
        Q_strict_avg = np.mean(Q[strict_mask])
        return np.where(Q <= Q_strict_avg)[0]

    strict_mask = _strict_baseflow(Q_filled, quantile)

    # Need at least some strict days to compute threshold
    if strict_mask.sum() == 0:
        result = np.zeros(n, dtype=bool)
        result[nan_mask] = False
        return result

    valid_indices = _flow_thresh(Q_filled, strict_mask)
    strict_indy = np.where(strict_mask)[0]
    combined_indy = np.union1d(valid_indices, strict_indy)

    result = np.zeros(n, dtype=bool)
    result[combined_indy] = True
    result[nan_mask] = False  # NaN days are never baseflow days
    return result


def separation_skill(hydrograph, SBT, basin_char, gw_hyd, flow, quantile=0.8):
    """Evaluate how well pybfs baseflow separation matches the modified strict baseflow reference.

    On days identified as baseflow-dominated by the modified strict method, observed
    streamflow Q is treated as the reference baseflow. Residuals Q - BF_bfs measure
    how well pybfs matches this reference.

    Parameters
    ----------
    hydrograph : pd.DataFrame
        Observed streamflow with columns 'Date' (datetime) and 'Streamflow' (m³/day).
    SBT : pd.DataFrame
        Baseflow table from base_table() with columns ['Xb','Z','S','Q'].
    basin_char : list
        [area, lb, x1, wb, por] — same as bfs().
    gw_hyd : list
        [alpha, beta, ks, kb, kz] — same as bfs().
    flow : list
        [qthresh, rs, rb1, rb2, prec, fr4rise] — same as bfs().
    quantile : float, optional
        Quantile value used in modified strict baseflow (default 0.8).

    Returns
    -------
    skill_df : pd.DataFrame
        Daily record with columns: Date, Q, BF_strict, BF_bfs, RES.
    metrics : dict
        Overall metrics: RMSE, MAE, n_days, frac_strict.
    """
    # Step 1: Compute modified strict baseflow mask
    bf_strict = modified_strict_baseflow(hydrograph, quantile=quantile)

    # Step 2: Run BFS separation
    bfs_result = bfs(hydrograph, SBT, basin_char, gw_hyd, flow)

    # Step 3: Build skill_df
    skill_df = pd.DataFrame({
        'Date': pd.to_datetime(hydrograph['Date'].values),
        'Q': np.array(hydrograph['Streamflow'], dtype=float),
        'BF_strict': bf_strict,
        'BF_bfs': bfs_result['Baseflow'].values,
        'RES': np.nan,
    })

    # Step 4: Compute residuals on strict-baseflow days
    mask = skill_df['BF_strict'] & skill_df['Q'].notna() & skill_df['BF_bfs'].notna()
    skill_df.loc[mask, 'RES'] = skill_df.loc[mask, 'Q'] - skill_df.loc[mask, 'BF_bfs']

    # Step 5: Compute overall metrics
    res_valid = skill_df.loc[mask, 'RES'].values
    n_days = int(mask.sum())
    total_days = len(skill_df)
    if n_days > 0:
        rmse = float(np.sqrt(np.mean(res_valid ** 2)))
        mae = float(np.mean(np.abs(res_valid)))
    else:
        rmse = np.nan
        mae = np.nan

    metrics = {
        'RMSE': rmse,
        'MAE': mae,
        'n_days': n_days,
        'frac_strict': n_days / total_days if total_days > 0 else np.nan,
    }

    return skill_df, metrics


def forecast_skill(hydrograph, SBT, basin_char, gw_hyd, flow,
                   min_days=30, max_days=90, min_sat=0.9, train_days=365,
                   quantile=0.8):
    """Evaluate pybfs forecast skill during naturally-occurring baseflow-dominated periods.

    Identifies contiguous sequences of baseflow-dominated days (via modified_strict_baseflow),
    runs a pybfs forecast initialized from observed conditions before each sequence, and
    computes error metrics comparing the forecasted baseflow to observed streamflow.

    Parameters
    ----------
    hydrograph : pd.DataFrame
        Observed streamflow with columns 'Date' (datetime) and 'Streamflow' (m³/day).
    SBT : pd.DataFrame
        Baseflow table from base_table() with columns ['Xb','Z','S','Q'].
    basin_char : list
        [area, lb, x1, wb, por] — same as bfs().
    gw_hyd : list
        [alpha, beta, ks, kb, kz] — same as bfs().
    flow : list
        [qthresh, rs, rb1, rb2, prec, fr4rise] — same as bfs().
    min_days : int, optional
        Minimum sequence length in days to qualify (default 30).
    max_days : int, optional
        Maximum forecast period — sequences longer than this are truncated (default 90).
    min_sat : float, optional
        Minimum fraction of baseflow days within a sequence to qualify (default 0.9).
    train_days : int, optional
        Number of days of observed data before the sequence used to initialize the
        forecast. Uses all available data if fewer days exist (default 365).
    quantile : float, optional
        Quantile value used in modified strict baseflow (default 0.8).

    Returns
    -------
    skill_df : pd.DataFrame
        Daily record with columns: Date, Q, BF, SEQ, FC, RES.
    summary_df : pd.DataFrame
        One row per qualifying sequence with columns: SEQ, LEN, SAT, RMSE, MAE.
    metrics : dict
        Overall metrics across all sequences: overall_RMSE, overall_MAE.
    """
    hydrograph = hydrograph.reset_index(drop=True)
    dates = pd.to_datetime(hydrograph['Date'].values)
    Q = np.array(hydrograph['Streamflow'], dtype=float)
    n = len(hydrograph)

    # Step 1: Compute modified strict baseflow mask
    bf_mask = modified_strict_baseflow(hydrograph, quantile=quantile)

    # Step 2: Initialize skill_df
    skill_df = pd.DataFrame({
        'Date': dates,
        'Q': Q,
        'BF': bf_mask,
        'SEQ': np.nan,
        'FC': np.nan,
        'RES': np.nan,
    })

    # Step 3: Find contiguous runs of True in bf_mask
    sequences = []  # list of (start_idx, end_idx) inclusive
    i = 0
    while i < n:
        if bf_mask[i]:
            j = i
            while j < n and bf_mask[j]:
                j += 1
            run_len = j - i
            run_sat = run_len / run_len  # always 1.0 for contiguous True runs
            if run_len >= min_days and run_sat >= min_sat:
                sequences.append((i, j - 1))  # inclusive end
            i = j
        else:
            i += 1

    # Step 4: Forecast each qualifying sequence
    summary_rows = []
    seq_num = 0

    for (start_idx, end_idx) in sequences:
        seq_num += 1
        seq_len = min(end_idx - start_idx + 1, max_days)
        actual_end_idx = start_idx + seq_len - 1

        # Assign sequence number in skill_df
        skill_df.loc[start_idx:actual_end_idx, 'SEQ'] = seq_num

        # Determine training slice
        train_start_idx = max(0, start_idx - train_days)
        train_slice = hydrograph.iloc[train_start_idx:start_idx].reset_index(drop=True)

        if len(train_slice) == 0:
            # No prior data — skip forecasting this sequence
            skill_df.loc[start_idx:actual_end_idx, 'SEQ'] = np.nan
            seq_num -= 1
            continue

        # Run bfs on training slice to get initial conditions
        train_bfs = bfs(train_slice, SBT, basin_char, gw_hyd, flow)
        last = train_bfs.iloc[-1]
        ini = (
            last['X'],
            last['Zb.L'],
            last['Zs.L'],
            last['StBase'],
            last['StSur'],
            last['SurfaceFlow'],
            last['Baseflow'],
            last['Rech'],
        )

        # Build forecast input DataFrame (lowercase columns as expected by forecast())
        forecast_dates = dates[start_idx: start_idx + seq_len]
        forecast_input = pd.DataFrame({
            'date': forecast_dates,
            'streamflow': np.nan,
        })

        # Run forecast
        fc_result = forecast(forecast_input, SBT, basin_char, gw_hyd, flow, ini)

        # Populate FC in skill_df
        fc_values = fc_result['Baseflow'].values
        skill_df.loc[start_idx: start_idx + seq_len - 1, 'FC'] = fc_values

        # Compute residuals
        skill_df.loc[start_idx: start_idx + seq_len - 1, 'RES'] = (
            skill_df.loc[start_idx: start_idx + seq_len - 1, 'Q'].values
            - fc_values
        )

        # Compute per-sequence metrics
        res_vals = skill_df.loc[start_idx: start_idx + seq_len - 1, 'RES'].dropna().values
        sat = float(bf_mask[start_idx: start_idx + seq_len].sum()) / seq_len
        if len(res_vals) > 0:
            seq_rmse = float(np.sqrt(np.mean(res_vals ** 2)))
            seq_mae = float(np.mean(np.abs(res_vals)))
        else:
            seq_rmse = np.nan
            seq_mae = np.nan

        summary_rows.append({
            'SEQ': seq_num,
            'LEN': seq_len,
            'SAT': sat,
            'RMSE': seq_rmse,
            'MAE': seq_mae,
        })

    # Step 5: Build summary_df
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
    else:
        summary_df = pd.DataFrame(columns=['SEQ', 'LEN', 'SAT', 'RMSE', 'MAE'])

    # Step 6: Overall metrics
    all_res = skill_df['RES'].dropna().values
    if len(all_res) > 0:
        overall_rmse = float(np.sqrt(np.mean(all_res ** 2)))
        overall_mae = float(np.mean(np.abs(all_res)))
    else:
        overall_rmse = np.nan
        overall_mae = np.nan

    metrics = {
        'overall_RMSE': overall_rmse,
        'overall_MAE': overall_mae,
    }

    return skill_df, summary_df, metrics
