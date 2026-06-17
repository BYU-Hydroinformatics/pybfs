#!/usr/bin/env python3
"""
Test script for skill assessment functions in pybfs.

Runs all three skill functions on site 12167000 and prints results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pybfs

# ── Configuration ─────────────────────────────────────────────────────────────
site_number = 12167000
quantile = 0.8  # Quantile for modified_strict_baseflow (e.g., 0.25 for 25th percentile)

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
streamflow_data = pd.read_csv('bfs/12167000.csv')
streamflow_data = streamflow_data.rename(columns={'mean_daily_streamflow': 'Streamflow'})
streamflow_data['Date'] = pd.to_datetime(streamflow_data['Date'])
n_zero = (streamflow_data['Streamflow'] == 0).sum()
if n_zero:
    print(f"Masking {n_zero} zero-flow days as NaN")
    streamflow_data.loc[streamflow_data['Streamflow'] == 0, 'Streamflow'] = np.nan

bfs_params = pd.read_csv('bfs_params_12167000_python.csv')
bfs_params = bfs_params.rename(columns={'tmp.site': 'site_no', 'tmp.area': 'AREA'})
basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params, site_number)

lb, x1, wb, por = basin_char[1], basin_char[2], basin_char[3], basin_char[4]
alpha, beta, ks, kb, kz = gw_hyd
SBT = pybfs.base_table(lb, x1, wb, beta, kb, streamflow_data, por)

print(f"Record length: {len(streamflow_data)} days "
      f"({streamflow_data['Date'].iloc[0].date()} to {streamflow_data['Date'].iloc[-1].date()})\n")

# ── Phase 1: modified_strict_baseflow ─────────────────────────────────────────
print("=" * 60)
print("PHASE 1: modified_strict_baseflow")
print("=" * 60)

bf_mask = pybfs.modified_strict_baseflow(streamflow_data, quantile=quantile)

total_days = len(bf_mask)
n_strict = bf_mask.sum()
print(f"Total days:          {total_days}")
print(f"Strict-baseflow days: {n_strict} ({n_strict / total_days:.1%})")
print(f"Mean Q on all days:       {streamflow_data['Streamflow'].mean():.2f} m³/day")
print(f"Mean Q on strict days:    {streamflow_data['Streamflow'][bf_mask].mean():.2f} m³/day")

fig, ax = plt.subplots(figsize=(11, 4))
ax.plot(streamflow_data['Date'], streamflow_data['Streamflow'],
        color='steelblue', lw=0.8, label='Streamflow')
ax.scatter(streamflow_data['Date'][bf_mask], streamflow_data['Streamflow'][bf_mask],
           s=6, color='orangered', label='Strict baseflow days')
ax.set_yscale('log')
ax.set_ylabel('Q (m³/day)')
ax.set_title(f'Site {site_number} — modified_strict_baseflow')
ax.legend(loc='best')
fig.tight_layout()
plt.show()

# ── Phase 2: separation_skill ─────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PHASE 2: separation_skill")
print("=" * 60)

skill_df, metrics = pybfs.separation_skill(
    streamflow_data, SBT, basin_char, gw_hyd, flow, quantile=quantile
)

mask2 = skill_df['BF_strict'].astype(bool) & skill_df['RES'].notna()
mean_Q2 = skill_df.loc[mask2, 'Q'].mean()
sum_Q2 = skill_df.loc[mask2, 'Q'].sum()
nrmse2 = 100 * metrics['RMSE'] / mean_Q2 if mean_Q2 > 0 else np.nan
pbias2 = 100 * skill_df.loc[mask2, 'RES'].sum() / sum_Q2 if sum_Q2 > 0 else np.nan

print(f"RMSE:         {metrics['RMSE']:.2f} m³/day")
print(f"MAE:          {metrics['MAE']:.2f} m³/day")
print(f"NRMSE:        {nrmse2:.1f}%  (RMSE / mean Q on strict days)")
print(f"PBIAS:        {pbias2:+.1f}%  (positive = pybfs under-predicts)")
print(f"Strict days:  {metrics['n_days']} ({metrics['frac_strict']:.1%} of record)")
print("\nskill_df head:")
print(skill_df.head(10).to_string(index=False))

fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax2a.plot(skill_df['Date'], skill_df['Q'], color='steelblue', lw=0.8, label='Q (observed)')
ax2a.plot(skill_df['Date'], skill_df['BF_bfs'], color='darkgreen', lw=1.0, label='BF_bfs')
strict = skill_df['BF_strict'].astype(bool)
ax2a.scatter(skill_df['Date'][strict], skill_df['Q'][strict],
             s=6, color='orangered', label='Strict baseflow days')
ax2a.set_yscale('log')
ax2a.set_ylabel('Q (m³/day)')
ax2a.set_title(f'Site {site_number} — separation_skill '
               f'(RMSE={metrics["RMSE"]:.0f}, MAE={metrics["MAE"]:.0f}, '
               f'NRMSE={nrmse2:.1f}%, PBIAS={pbias2:+.1f}%)')
ax2a.legend(loc='best')

ax2b.bar(skill_df['Date'], skill_df['RES'], width=1.0, color='purple')
ax2b.axhline(0, color='k', lw=0.6)
ax2b.set_ylabel('Residual (Q − BF_bfs)')
ax2b.set_xlabel('Date')
fig2.tight_layout()
plt.show()

# ── Phase 3: forecast_skill ───────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PHASE 3: forecast_skill")
print("=" * 60)

skill_df_fc, summary_df, fc_metrics = pybfs.forecast_skill(
    streamflow_data, SBT, basin_char, gw_hyd, flow,
    min_days=30, max_days=90, min_sat=0.9, train_days=365, quantile=quantile
)

valid_fc = skill_df_fc['FC'].notna() & skill_df_fc['Q'].notna()
mean_Q_fc = skill_df_fc.loc[valid_fc, 'Q'].mean()
sum_Q_fc = skill_df_fc.loc[valid_fc, 'Q'].sum()
overall_nrmse = 100 * fc_metrics['overall_RMSE'] / mean_Q_fc if mean_Q_fc > 0 else np.nan
overall_pbias = (100 * skill_df_fc.loc[valid_fc, 'RES'].sum() / sum_Q_fc
                 if sum_Q_fc > 0 else np.nan)

nrmse_seq, pbias_seq = [], []
for s in summary_df['SEQ']:
    sub = skill_df_fc[(skill_df_fc['SEQ'] == s) & skill_df_fc['FC'].notna()]
    mq, sq = sub['Q'].mean(), sub['Q'].sum()
    rmse_s = summary_df.loc[summary_df['SEQ'] == s, 'RMSE'].iloc[0]
    nrmse_seq.append(100 * rmse_s / mq if mq > 0 else np.nan)
    pbias_seq.append(100 * sub['RES'].sum() / sq if sq > 0 else np.nan)
summary_df['NRMSE%'] = nrmse_seq
summary_df['PBIAS%'] = pbias_seq

print(f"Number of qualifying sequences: {len(summary_df)}")
print(f"Overall RMSE:  {fc_metrics['overall_RMSE']:.2f} m³/day")
print(f"Overall MAE:   {fc_metrics['overall_MAE']:.2f} m³/day")
print(f"Overall NRMSE: {overall_nrmse:.1f}%")
print(f"Overall PBIAS: {overall_pbias:+.1f}%")
print("\nSequence summary:")
print(summary_df.to_string(index=False))

if len(summary_df) > 0:
    seq1_days = skill_df_fc[skill_df_fc['SEQ'] == 1]
    print(f"\nSequence 1: {len(seq1_days)} days, "
          f"{seq1_days['Date'].iloc[0].date()} to {seq1_days['Date'].iloc[-1].date()}")
    print(seq1_days[['Date', 'Q', 'FC', 'RES']].head(10).to_string(index=False))

    fig_ov, ax_ov = plt.subplots(figsize=(11, 4))
    ax_ov.plot(skill_df_fc['Date'], skill_df_fc['Q'],
               color='steelblue', lw=0.8, label='Q (observed)')
    bf_ov = skill_df_fc['BF'].astype(bool)
    ax_ov.scatter(skill_df_fc['Date'][bf_ov], skill_df_fc['Q'][bf_ov],
                  s=6, color='orangered', label='Strict baseflow days', zorder=3)
    ax_ov.set_yscale('log')
    ax_ov.set_ylabel('Q (m³/day)')
    ax_ov.set_xlabel('Date')
    ax_ov.set_title(f'Site {site_number} — forecast_skill sequences '
                    f'(overall NRMSE={overall_nrmse:.1f}%, '
                    f'PBIAS={overall_pbias:+.1f}%)')
    for seq in sorted(summary_df['SEQ'].unique()):
        sub = skill_df_fc[skill_df_fc['SEQ'] == seq]
        t0, t1 = sub['Date'].iloc[0], sub['Date'].iloc[-1]
        ax_ov.axvspan(t0, t1, color='gray', alpha=0.3)
        y_top = ax_ov.get_ylim()[1]
        ax_ov.text(t0 + (t1 - t0) / 2, y_top, f'Seq {int(seq)}',
                   ha='center', va='top', fontsize=9, color='black')
    ax_ov.legend(loc='lower left')
    fig_ov.tight_layout()
    plt.show()

    n_seq = int(summary_df['SEQ'].max())
    fig3, axes = plt.subplots(n_seq, 1, figsize=(10, 4.2 * n_seq),
                              squeeze=False, constrained_layout=True)
    for i, seq in enumerate(sorted(summary_df['SEQ'].unique())):
        ax = axes[i, 0]
        sub = skill_df_fc[skill_df_fc['SEQ'] == seq]
        ax.plot(sub['Date'], sub['Q'], color='steelblue', lw=1.0, label='Q (observed)')
        ax.plot(sub['Date'], sub['FC'], color='crimson', lw=1.0, label='Forecast')
        row = summary_df[summary_df['SEQ'] == seq].iloc[0]
        ax.set_title(f'Sequence {int(seq)} — LEN={int(row["LEN"])}, '
                     f'RMSE={row["RMSE"]:.0f}, MAE={row["MAE"]:.0f}, '
                     f'NRMSE={row["NRMSE%"]:.1f}%, PBIAS={row["PBIAS%"]:+.1f}%')
        ax.set_ylabel('Q (m³/day)')
        ax.legend(loc='best')
    axes[-1, 0].set_xlabel('Date')
    for ax in axes[:, 0]:
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(30)
            lbl.set_ha('right')
    fig3.suptitle(f'Site {site_number} — forecast_skill '
                  f'(overall RMSE={fc_metrics["overall_RMSE"]:.0f}, '
                  f'NRMSE={overall_nrmse:.1f}%, PBIAS={overall_pbias:+.1f}%)')
    plt.show()

print("\n=== Done ===")
