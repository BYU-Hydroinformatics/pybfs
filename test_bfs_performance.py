#!/usr/bin/env python3
"""Test script to compare BFS performance with and without NUMBA JIT"""

import time
import pandas as pd
import numpy as np
import pybfs

def test_bfs_performance():
    """Compare performance of BFS with and without JIT"""
    
    # Load test data
    print("Loading test data...")
    streamflow_data = pd.read_csv('docs/files/2312200_data.csv')
    bfs_params_usgs = pd.read_csv('docs/files/bfs_params_50.csv')
    
    site_number = 2312200
    basin_char, gw_hyd, flow = pybfs.get_values_for_site(bfs_params_usgs, site_number)
    
    area, lb, x1, wb, por = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4]
    alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]
    
    print("Generating baseflow table...")
    SBT = pybfs.base_table(lb, x1, wb, beta, kb, streamflow_data, por)
    
    # Test without JIT
    print("\n" + "="*60)
    print("Testing BFS WITHOUT JIT compilation...")
    print("="*60)
    start_time = time.time()
    try:
        result_no_jit = pybfs.bfs(streamflow_data, SBT, basin_char, gw_hyd, flow, use_jit=False)
        time_no_jit = time.time() - start_time
        print(f"Time without JIT: {time_no_jit:.4f} seconds")
    except TypeError:
        # If use_jit parameter doesn't exist yet, just run normally
        result_no_jit = pybfs.bfs(streamflow_data, SBT, basin_char, gw_hyd, flow)
        time_no_jit = time.time() - start_time
        print(f"Time (original implementation): {time_no_jit:.4f} seconds")
        print("Note: NUMBA JIT not yet implemented - this is baseline timing")
    
    # Test with JIT (if implemented)
    print("\n" + "="*60)
    print("Testing BFS WITH JIT compilation...")
    print("="*60)
    try:
        # Warm-up run
        print("Warming up JIT compilation (first run may be slower)...")
        start_time = time.time()
        result_jit_warmup = pybfs.bfs(streamflow_data, SBT, basin_char, gw_hyd, flow, use_jit=True)
        time_jit_warmup = time.time() - start_time
        print(f"Time with JIT (warm-up): {time_jit_warmup:.4f} seconds")
        
        # Actual timing run
        print("Running JIT-compiled version (after warm-up)...")
        start_time = time.time()
        result_jit = pybfs.bfs(streamflow_data, SBT, basin_char, gw_hyd, flow, use_jit=True)
        time_jit = time.time() - start_time
        print(f"Time with JIT: {time_jit:.4f} seconds")
        
        # Compare results
        print("\n" + "="*60)
        print("Results Comparison:")
        print("="*60)
        speedup = time_no_jit / time_jit
        print(f"Speedup: {speedup:.2f}x faster with JIT")
        print(f"Time saved: {time_no_jit - time_jit:.4f} seconds")
        
        # Verify results are similar
        print("\n" + "="*60)
        print("Verifying results match...")
        print("="*60)
        
        # Compare key columns
        columns_to_compare = ['Qsim', 'Baseflow', 'SurfaceFlow', 'DirectRunoff']
        max_diff = 0.0
        for col in columns_to_compare:
            diff = np.nanmax(np.abs(result_no_jit[col] - result_jit[col]))
            max_diff = max(max_diff, diff)
            print(f"Max difference in {col}: {diff:.6e}")
        
        if max_diff < 1e-6:
            print("\n✓ Results match within numerical precision!")
        else:
            print(f"\n⚠ Warning: Results differ by up to {max_diff:.6e}")
        
        return {
            'time_no_jit': time_no_jit,
            'time_jit_warmup': time_jit_warmup,
            'time_jit': time_jit,
            'speedup': speedup,
            'max_diff': max_diff
        }
        
    except TypeError:
        print("NUMBA JIT not yet implemented in bfs() function")
        print("Please implement NUMBA JIT first (see previous suggestions)")
        return {
            'time_no_jit': time_no_jit,
            'time_jit': None,
            'speedup': None
        }

if __name__ == '__main__':
    results = test_bfs_performance()
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    print(f"Original time: {results['time_no_jit']:.4f}s")
    if results['time_jit'] is not None:
        print(f"JIT time (after warm-up): {results['time_jit']:.4f}s")
        print(f"Speedup: {results['speedup']:.2f}x")
    else:
        print("JIT not available - implement NUMBA JIT first")

