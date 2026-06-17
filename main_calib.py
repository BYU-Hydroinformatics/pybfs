#!/usr/bin/env python3
"""
Calibration script for site 12167000 - Python vs R comparison

This script:
1. Loads streamflow data for gage 12167000
2. Runs the Python calibration process
3. Compares results with R calibration parameters
"""

import pandas as pd
import numpy as np
import pybfs


def main():
    """Main execution function for calibration testing"""

    # Site information
    site_id = "12167000"
    site_area = 6.71e+08  # m² (from siteinfo file)

    # Load streamflow data
    print("Loading streamflow data...")
    streamflow_data = pd.read_csv('bfs/12167000.csv', encoding='utf-8-sig')
    streamflow_data['Date'] = pd.to_datetime(streamflow_data['Date'], format='%m/%d/%Y')
    
    # Extract streamflow and dates
    tmp_q = streamflow_data['mean_daily_streamflow'].values
    dys = streamflow_data['Date'].values

    print(f"Loaded {len(tmp_q)} days of streamflow data")
    print(f"Date range: {dys[0]} to {dys[-1]}")
    print(f"Mean streamflow: {np.nanmean(tmp_q):.2f} m³/day")

    # Load R calibration parameters for comparison
    print("\nLoading R calibration parameters...")
    r_params = pd.read_csv('bfs/out/site_sum/bfs_params_12167000.csv')
    r_bff = pd.read_csv('bfs/out/site_sum/bff_12167000.csv')

    print(f"\nR calibration parameters for site {site_id}:")
    print(f"  AREA: {r_params['tmp.area'].iloc[0]}")
    print(f"  Lb: {r_params['Lb'].iloc[0]:.6f}")
    print(f"  X1: {r_params['X1'].iloc[0]:.6f}")
    print(f"  Wb: {r_params['Wb'].iloc[0]:.6f}")
    print(f"  POR: {r_params['POR'].iloc[0]:.6f}")
    print(f"  ALPHA: {r_params['ALPHA'].iloc[0]:.6f}")
    print(f"  BETA: {r_params['BETA'].iloc[0]:.6f}")
    print(f"  Ks: {r_params['Ks'].iloc[0]:.6f}")
    print(f"  Kb: {r_params['Kb'].iloc[0]:.6f}")
    print(f"  Kz: {r_params['Kz'].iloc[0]:.6f}")
    print(f"  Qthresh: {r_params['Qthresh'].iloc[0]:.6f}")
    print(f"  Rs: {r_params['Rs'].iloc[0]:.6f}")
    print(f"  Rb1: {r_params['Rb1'].iloc[0]:.6f}")
    print(f"  Rb2: {r_params['Rb2'].iloc[0]:.6f}")
    print(f"  Prec: {r_params['Prec'].iloc[0]:.6f}")
    print(f"  Frac4Rise: {r_params['Frac4Rise'].iloc[0]:.6f}")
    print(f"  Error: {r_params['Error'].iloc[0]:.6f}")
    print(f"  BFF: {r_params['BFF'].iloc[0]:.6f}")

    print(f"\nR flow fractions:")
    print(f"  Qmean: {r_bff['Qmean'].iloc[0]:.6f}")
    print(f"  BFF: {r_bff['BFF'].iloc[0]:.6f}")
    print(f"  SFF: {r_bff['SFF'].iloc[0]:.6f}")
    print(f"  DRF: {r_bff['DRF'].iloc[0]:.6f}")
    print(f"  Error: {r_bff['Error'].iloc[0]:.6f}")

    # Run Python calibration
    print("\n" + "="*60)
    print("Running Python calibration...")
    print("="*60)
    print("This may take several minutes...")

    try:
        bf_params, bff, ci_table, bfs_out = pybfs.bfs_calibrate(
            tmp_site=site_id,
            tmp_area=site_area,
            tmp_q=tmp_q,
            dys=dys
        )

        if bf_params is None:
            print("\nERROR: Calibration failed - flow_metrics returned invalid values")
            return

        print("\n" + "="*60)
        print("PYTHON CALIBRATION COMPLETE")
        print("="*60)

        # Display calibrated parameters
        print("\nPython calibrated parameters:")
        print(f"  AREA: {bf_params['tmp.area'].iloc[0]}")
        print(f"  Lb: {bf_params['Lb'].iloc[0]:.6f}")
        print(f"  X1: {bf_params['X1'].iloc[0]:.6f}")
        print(f"  Wb: {bf_params['Wb'].iloc[0]:.6f}")
        print(f"  POR: {bf_params['POR'].iloc[0]:.6f}")
        print(f"  ALPHA: {bf_params['ALPHA'].iloc[0]:.6f}")
        print(f"  BETA: {bf_params['BETA'].iloc[0]:.6f}")
        print(f"  Ks: {bf_params['Ks'].iloc[0]:.6f}")
        print(f"  Kb: {bf_params['Kb'].iloc[0]:.6f}")
        print(f"  Kz: {bf_params['Kz'].iloc[0]:.6f}")
        print(f"  Qthresh: {bf_params['Qthresh'].iloc[0]:.6f}")
        print(f"  Rs: {bf_params['Rs'].iloc[0]:.6f}")
        print(f"  Rb1: {bf_params['Rb1'].iloc[0]:.6f}")
        print(f"  Rb2: {bf_params['Rb2'].iloc[0]:.6f}")
        print(f"  Prec: {bf_params['Prec'].iloc[0]:.6f}")
        print(f"  Frac4Rise: {bf_params['Frac4Rise'].iloc[0]:.6f}")
        print(f"  Error: {bf_params['Error'].iloc[0]:.6f}")
        print(f"  BFF: {bf_params['BFF'].iloc[0]:.6f}")

        # Display flow fractions
        print("\nPython flow fractions:")
        print(f"  Qmean: {bff['Qmean'].iloc[0]:.6f}")
        print(f"  BFF (Baseflow Fraction): {bff['BFF'].iloc[0]:.6f}")
        print(f"  SFF (Surface Flow Fraction): {bff['SFF'].iloc[0]:.6f}")
        print(f"  DRF (Direct Runoff Fraction): {bff['DRF'].iloc[0]:.6f}")
        print(f"  Error: {bff['Error'].iloc[0]:.6f}")

        # Compare with R results
        print("\n" + "="*60)
        print("COMPARISON: PYTHON vs R CALIBRATION PARAMETERS")
        print("="*60)

        param_names = ['Lb', 'X1', 'Wb', 'POR', 'ALPHA', 'BETA', 'Ks', 'Kb', 'Kz',
                      'Qthresh', 'Rs', 'Rb1', 'Rb2', 'Prec', 'Frac4Rise', 'Error', 'BFF']

        print(f"\n{'Parameter':<15} {'R Value':<15} {'Python Value':<15} {'Difference':<15} {'% Diff':<15}")
        print("-" * 75)

        for param in param_names:
            r_val = r_params[param].iloc[0]
            py_val = bf_params[param].iloc[0]
            diff = py_val - r_val
            pct_diff = (diff / r_val * 100) if r_val != 0 else 0
            print(f"{param:<15} {r_val:<15.6f} {py_val:<15.6f} {diff:<15.6f} {pct_diff:<15.2f}%")

        # Summary
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        
        # Calculate relative differences for key parameters
        key_params = ['Lb', 'X1', 'Wb', 'ALPHA', 'BETA', 'Ks', 'Kb', 'Kz']
        max_diff = 0
        max_param = None
        
        for param in key_params:
            r_val = r_params[param].iloc[0]
            py_val = bf_params[param].iloc[0]
            if r_val != 0:
                pct_diff = abs((py_val - r_val) / r_val * 100)
                if pct_diff > max_diff:
                    max_diff = pct_diff
                    max_param = param

        print(f"\nMaximum relative difference in key parameters: {max_diff:.2f}% ({max_param})")
        
        if max_diff < 1:
            print("✓ Python and R calibration results are very close!")
        elif max_diff < 5:
            print("⚠ Python and R calibration results are reasonably close")
        elif max_diff < 10:
            print("⚠ Python and R calibration results show some differences")
        else:
            print("✗ Python and R calibration results differ significantly")
            print("  This may be due to:")
            print("  - Different optimization algorithms or convergence criteria")
            print("  - Different random initialization")
            print("  - Numerical precision differences between R and Python")
            print("  - Different optimization iteration limits")

        # Compare flow fractions
        print("\n" + "="*60)
        print("FLOW FRACTIONS COMPARISON")
        print("="*60)
        print(f"\n{'Metric':<15} {'R Value':<15} {'Python Value':<15} {'Difference':<15} {'% Diff':<15}")
        print("-" * 75)
        
        for metric in ['Qmean', 'BFF', 'SFF', 'DRF', 'Error']:
            r_val = r_bff[metric].iloc[0]
            py_val = bff[metric].iloc[0]
            diff = py_val - r_val
            pct_diff = (diff / r_val * 100) if r_val != 0 else 0
            print(f"{metric:<15} {r_val:<15.6f} {py_val:<15.6f} {diff:<15.6f} {pct_diff:<15.2f}%")

        # Save Python calibrated parameters in same format as R
        print("\n" + "="*60)
        print("SAVING PYTHON CALIBRATION RESULTS")
        print("="*60)
        
        # Save bfs_params (same format as R)
        bf_params_output = bf_params.copy()
        bf_params_output.to_csv(
            f'bfs_params_{site_id}_python.csv',
            index=False,
            float_format='%.6g'
        )
        print(f"Saved: bfs_params_{site_id}_python.csv")
        
        # Save bff (same format as R)
        bff_output = bff.copy()
        bff_output.to_csv(
            f'bff_{site_id}_python.csv',
            index=False,
            float_format='%.6g'
        )
        print(f"Saved: bff_{site_id}_python.csv")

        # Compare BFS results with R vs Python parameters
        print("\n" + "="*60)
        print("COMPARING BFS RESULTS: R vs PYTHON CALIBRATION PARAMETERS")
        print("="*60)
        
        try:
            import matplotlib.pyplot as plt
            
            # Prepare streamflow DataFrame
            streamflow_df = pd.DataFrame({
                'Date': pd.to_datetime(dys),
                'Streamflow': tmp_q
            })
            
            # Get flow metrics (already calculated during calibration)
            flow_metrics_result = pybfs.flow_metrics(tmp_q, timestep='day', fr4rise=0.05)
            flow = [
                flow_metrics_result[0],  # Qthresh
                flow_metrics_result[1],  # Rs
                flow_metrics_result[2],  # Rb1
                flow_metrics_result[3],  # Rb2
                flow_metrics_result[4],  # Prec
                flow_metrics_result[5]   # Fr4Rise
            ]
            
            # PYTHON CALIBRATED PARAMETERS
            print("\nRunning BFS with Python calibrated parameters...")
            py_basin_char = [
                site_area,
                bf_params['Lb'].iloc[0],
                bf_params['X1'].iloc[0],
                bf_params['Wb'].iloc[0],
                bf_params['POR'].iloc[0]
            ]
            py_gw_hyd = [
                bf_params['ALPHA'].iloc[0],
                bf_params['BETA'].iloc[0],
                bf_params['Ks'].iloc[0],
                bf_params['Kb'].iloc[0],
                bf_params['Kz'].iloc[0]
            ]
            py_SBT = pybfs.base_table(
                py_basin_char[1], py_basin_char[2], py_basin_char[3],
                py_gw_hyd[1], py_gw_hyd[3], streamflow_df, py_basin_char[4]
            )
            py_bfs_out = pybfs.bfs(streamflow_df, py_SBT, py_basin_char, py_gw_hyd, flow)
            
            # R CALIBRATED PARAMETERS
            print("Running BFS with R calibrated parameters...")
            r_basin_char = [
                site_area,
                r_params['Lb'].iloc[0],
                r_params['X1'].iloc[0],
                r_params['Wb'].iloc[0],
                r_params['POR'].iloc[0]
            ]
            r_gw_hyd = [
                r_params['ALPHA'].iloc[0],
                r_params['BETA'].iloc[0],
                r_params['Ks'].iloc[0],
                r_params['Kb'].iloc[0],
                r_params['Kz'].iloc[0]
            ]
            r_SBT = pybfs.base_table(
                r_basin_char[1], r_basin_char[2], r_basin_char[3],
                r_gw_hyd[1], r_gw_hyd[3], streamflow_df, r_basin_char[4]
            )
            r_bfs_out = pybfs.bfs(streamflow_df, r_SBT, r_basin_char, r_gw_hyd, flow)
            
            # Create comparison plot
            print("Creating comparison plot...")
            fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
            
            # Determine column names (handle both formats)
            qob_col = 'Qob' if 'Qob' in py_bfs_out.columns else 'Qob.L3'
            baseflow_col = 'Baseflow' if 'Baseflow' in py_bfs_out.columns else 'Baseflow.L3'
            surface_col = 'SurfaceFlow' if 'SurfaceFlow' in py_bfs_out.columns else 'SurfaceFlow.L3'
            direct_col = 'DirectRunoff' if 'DirectRunoff' in py_bfs_out.columns else 'DirectRunoff.L3'
            
            dates = py_bfs_out['Date'].values
            
            # Calculate total simulated flow for statistics (not plotted)
            py_qsim = py_bfs_out[baseflow_col].values + py_bfs_out[surface_col].values + py_bfs_out[direct_col].values
            r_qsim = r_bfs_out[baseflow_col].values + r_bfs_out[surface_col].values + r_bfs_out[direct_col].values
            
            # Plot 1: Observed Streamflow with R and Python Baseflow
            axes[0].plot(dates, py_bfs_out[qob_col].values, 'k-', linewidth=1.5, label='Observed Streamflow', alpha=0.7)
            axes[0].plot(dates, py_bfs_out[baseflow_col].values, 'b-', linewidth=1.5, 
                        label='Python Baseflow', alpha=0.7)
            axes[0].plot(dates, r_bfs_out[baseflow_col].values, 'r--', linewidth=1.5, 
                        label='R Baseflow', alpha=0.7)
            axes[0].set_ylabel('Flow\n(m³/day)', fontsize=10)
            axes[0].set_title('Observed Streamflow and Baseflow Comparison', fontsize=12, fontweight='bold')
            axes[0].grid(True, alpha=0.3)
            axes[0].legend()
            
            # Plot 2: Baseflow Comparison
            axes[1].plot(dates, py_bfs_out[baseflow_col].values, 'b-', linewidth=1.5, 
                        label='Python Baseflow', alpha=0.7)
            axes[1].plot(dates, r_bfs_out[baseflow_col].values, 'r--', linewidth=1.5, 
                        label='R Baseflow', alpha=0.7)
            axes[1].set_ylabel('Baseflow\n(m³/day)', fontsize=10)
            axes[1].set_title('Baseflow Comparison: Python vs R', fontsize=12, fontweight='bold')
            axes[1].grid(True, alpha=0.3)
            axes[1].legend()
            
            # Plot 3: Surface Flow Comparison
            axes[2].plot(dates, py_bfs_out[surface_col].values, 'b-', linewidth=1.5, 
                        label='Python Surface Flow', alpha=0.7)
            axes[2].plot(dates, r_bfs_out[surface_col].values, 'r--', linewidth=1.5, 
                        label='R Surface Flow', alpha=0.7)
            axes[2].set_ylabel('Surface Flow\n(m³/day)', fontsize=10)
            axes[2].set_xlabel('Date', fontsize=10)
            axes[2].set_title('Surface Flow Comparison: Python vs R', fontsize=12, fontweight='bold')
            axes[2].grid(True, alpha=0.3)
            axes[2].legend()
            
            plt.tight_layout()
            plot_filename = f'bfs_comparison_{site_id}.png'
            plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
            print(f"\nComparison plot saved as '{plot_filename}'")
            
            # Print summary statistics
            print("\n" + "="*60)
            print("SUMMARY STATISTICS")
            print("="*60)
            
            py_bff_actual = np.nansum(py_bfs_out[baseflow_col]) / np.nansum(py_bfs_out[qob_col])
            r_bff_actual = np.nansum(r_bfs_out[baseflow_col]) / np.nansum(r_bfs_out[qob_col])
            
            print(f"\nBaseflow Fraction (BFF) from time series:")
            print(f"  Python: {py_bff_actual:.4f}")
            print(f"  R:      {r_bff_actual:.4f}")
            print(f"  Difference:  {py_bff_actual - r_bff_actual:.4f} ({((py_bff_actual - r_bff_actual) / r_bff_actual * 100):.2f}%)")
            
            # Calculate mean absolute error for simulated flow
            py_mae = np.nanmean(np.abs(py_qsim - py_bfs_out[qob_col].values))
            r_mae = np.nanmean(np.abs(r_qsim - r_bfs_out[qob_col].values))
            
            print(f"\nMean Absolute Error (Simulated vs Observed):")
            print(f"  Python: {py_mae:.2f} m³/day")
            print(f"  R:      {r_mae:.2f} m³/day")
            
            plt.show()
            
        except Exception as e:
            print(f"\nERROR during plotting: {e}")
            import traceback
            traceback.print_exc()

    except Exception as e:
        print(f"\nERROR during calibration: {e}")
        import traceback
        traceback.print_exc()
        return


if __name__ == "__main__":
    main()

