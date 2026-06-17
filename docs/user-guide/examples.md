# Usage Examples

This page demonstrates practical applications of PyBFS for baseflow separation and forecasting.

## Interactive notebook (Google Colab)

If you want to run these examples in a browser (no local setup), use the Colab notebook:

<a href="https://colab.research.google.com/github/BYU-Hydroinformatics/pybfs/blob/main/notebooks/pybfs_forecast.ipynb" target="_"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

At the first part of the notebook, you will be prompted to upload a streamflow file and a parameters file. You can download the required files here:

Streamflow: [2312200_data.csv](../files/2312200_data.csv)<br>
Parameters: [bfs_params_50.csv](../files/bfs_params_50.csv)

The notebook includes code for a basic baseflow separation analysis and a short-term baseflow forecast.

## Baseflow Separation Example

This example demonstrates how to use PyBFS to separate baseflow from streamflow data for a specific site. This assumes you have already installed PyBFS and have the necessary data files and are running within a Python environment.


```python
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pybfs as bfs

# Load streamflow data
streamflow_data = pd.read_csv('docs/files/2312200_data.csv')

# Load site parameters
bfs_params_usgs = pd.read_csv('docs/files/bfs_params_50.csv')

# Get parameters for specific site
site_number = 2312200
basin_char, gw_hyd, flow = bfs.get_values_for_site(bfs_params_usgs, site_number)

# Extract basin characteristics
area, lb, x1, wb, por = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4]
ws = wb / 2

# Extract groundwater hydraulic parameters
alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]

# Extract flow metrics
qthresh, rs, rb1, rb2, prec, fr4rise = flow[0], flow[1], flow[2], flow[3], flow[4], flow[5]

print(f"Basin characteristics:")
print(f"  Area: {area}, Length: {lb}, Width: {wb}")
print(f"  Porosity: {por}")

# Generate baseflow table
SBT = bfs.base_table(lb, x1, wb, beta, kb, streamflow_data, por)
print(f"Baseflow table generated with {len(SBT)} rows")

# Run PyBFS
result = bfs.PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
print(f"PyBFS completed for {len(result)} time steps")

# Display summary statistics
print("\n=== Results Summary ===")
print(f"Total observed flow: {result['Qob'].sum():.2f}")
print(f"Total simulated flow: {result['Qsim'].sum():.2f}")
print(f"Total baseflow: {result['Baseflow'].sum():.2f}")
print(f"Total surface flow: {result['SurfaceFlow'].sum():.2f}")
print(f"Total direct runoff: {result['DirectRunoff'].sum():.2f}")

# Plot results
bfs.plot_baseflow_simulation(streamflow_data, result)
```

## Forecasting Example

Once you have calibrated the model, you can create forecasts for future periods as shown below. Again, this assumes you have the necessary data files and are running within a Python environment.

```python
# Filter data for calibration period (Jan-Sep 2018)
start_date = '2018-01-01'
end_date = '2018-09-30'
streamflow_data['Date'] = pd.to_datetime(streamflow_data['Date'])
streamflow_data_filtered = streamflow_data[
    (streamflow_data['Date'] >= start_date) & (streamflow_data['Date'] <= end_date)
]

# Run PyBFS for calibration period
tmp2 = bfs.PyBFS(streamflow_data_filtered, SBT, basin_char, gw_hyd, flow)

# Extract initial conditions from last time step
Xi, Zbi, Zsi, StBi, StSi, Surflow, Baseflow, Rech = tmp2.iloc[-1][
    ['X', 'Zb.L', 'Zs.L', 'StBase', 'StSur', 'SurfaceFlow', 'Baseflow', 'Rech']
]
ini = (Xi, Zbi, Zsi, StBi, StSi, Surflow, Baseflow, Rech)

print(f"Initial conditions extracted from {tmp2.iloc[-1]['Date']}")

# Create forecast period (Oct-Nov 2018)
dates = pd.date_range(start="2018-10-01", end="2018-11-30", freq="D")
forecast_df = pd.DataFrame({
    "date": dates,
    "streamflow": np.nan
})

# Run forecast
f = bfs.forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
print(f"Forecast completed for {len(f)} time steps")

# Plot forecast
bfs.plot_forecast_baseflow(f)

# Plot forecast with observed data for comparison
forecast_start = '2018-10-01'
forecast_end = '2018-11-30'
streamflow_data_forecast = streamflow_data[
    (streamflow_data['Date'] >= forecast_start) & (streamflow_data['Date'] <= forecast_end)
]

bfs.plot_forecast_baseflow_streamflow(f, streamflow_data_forecast)
```

## Calibration Example

This example demonstrates how to calibrate PyBFS parameters automatically from streamflow data. Calibration is useful when you don't have pre-calibrated parameters or want to optimize parameters for a specific site. The calibration process follows the methodology described in the USGS BFS manual (see `refs/usgs_bfs_manual.pdf` for detailed information about the calibration algorithm and objective functions).

```python
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pybfs as bfs

# Load streamflow data
streamflow_data = pd.read_csv('your_streamflow_data.csv')
streamflow_data['Date'] = pd.to_datetime(streamflow_data['Date'])

# Site information
site_id = "12345678"
site_area = 1.5e8  # Drainage area in m²

# Extract streamflow and dates
tmp_q = streamflow_data['Streamflow'].values
dys = streamflow_data['Date'].values

print(f"Calibrating parameters for site {site_id}...")
print(f"Using {len(tmp_q)} days of streamflow data")
print(f"Date range: {dys[0]} to {dys[-1]}")

# Run calibration
# Note: This may take several minutes depending on data length
bf_params, bff, ci_table, bfs_out = bfs.bfs_calibrate(
    tmp_site=site_id,
    tmp_area=site_area,
    tmp_q=tmp_q,
    dys=dys
)

# Check if calibration succeeded
if bf_params is None:
    print("Calibration failed - check your data quality")
else:
    print("\n=== Calibration Complete ===")
    print(f"\nCalibrated Parameters:")
    print(f"  Lb (Basin Length): {bf_params['Lb'].iloc[0]:.2f} m")
    print(f"  X1 (Base Gradient): {bf_params['X1'].iloc[0]:.2f} m")
    print(f"  Wb (Basin Width): {bf_params['Wb'].iloc[0]:.2f} m")
    print(f"  POR (Porosity): {bf_params['POR'].iloc[0]:.4f}")
    print(f"  ALPHA (Surface Gradient): {bf_params['ALPHA'].iloc[0]:.6f}")
    print(f"  BETA (Base Exponent): {bf_params['BETA'].iloc[0]:.3f}")
    print(f"  Ks (Surface Conductivity): {bf_params['Ks'].iloc[0]:.6f} m/day")
    print(f"  Kb (Base Conductivity): {bf_params['Kb'].iloc[0]:.6f} m/day")
    print(f"  Kz (Vertical Conductivity): {bf_params['Kz'].iloc[0]:.6f} m/day")
    
    print(f"\nFlow Fractions:")
    print(f"  BFF (Baseflow Fraction): {bff['BFF'].iloc[0]:.3f}")
    print(f"  SFF (Surface Flow Fraction): {bff['SFF'].iloc[0]:.3f}")
    print(f"  DRF (Direct Runoff Fraction): {bff['DRF'].iloc[0]:.3f}")
    print(f"  Error: {bff['Error'].iloc[0]:.6f}")
    
    # Extract calibrated parameters for use in baseflow separation
    basin_char = [
        site_area,
        bf_params['Lb'].iloc[0],
        bf_params['X1'].iloc[0],
        bf_params['Wb'].iloc[0],
        bf_params['POR'].iloc[0]
    ]
    
    gw_hyd = [
        bf_params['ALPHA'].iloc[0],
        bf_params['BETA'].iloc[0],
        bf_params['Ks'].iloc[0],
        bf_params['Kb'].iloc[0],
        bf_params['Kz'].iloc[0]
    ]
    
    flow = [
        bf_params['Qthresh'].iloc[0],
        bf_params['Rs'].iloc[0],
        bf_params['Rb1'].iloc[0],
        bf_params['Rb2'].iloc[0],
        bf_params['Prec'].iloc[0],
        bf_params['Frac4Rise'].iloc[0]
    ]
    
    # Generate baseflow table with calibrated parameters
    SBT = bfs.base_table(basin_char[1], basin_char[2], basin_char[3],
                         gw_hyd[1], gw_hyd[3], streamflow_data, basin_char[4])
    
    # Run baseflow separation with calibrated parameters
    result = bfs.PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
    
    # Plot results
    bfs.plot_baseflow_simulation(streamflow_data, result)
    
    # Save calibrated parameters for future use
    bf_params.to_csv(f'bfs_params_{site_id}.csv', index=False)
    print(f"\nCalibrated parameters saved to bfs_params_{site_id}.csv")
```

### Using Calibrated Parameters

Once you have calibrated parameters, you can use them for baseflow separation on the same site or similar sites:

```python
# Load previously calibrated parameters
bf_params = pd.read_csv('bfs_params_12345678.csv')

# Extract parameters
basin_char = [
    bf_params['tmp.area'].iloc[0],
    bf_params['Lb'].iloc[0],
    bf_params['X1'].iloc[0],
    bf_params['Wb'].iloc[0],
    bf_params['POR'].iloc[0]
]

gw_hyd = [
    bf_params['ALPHA'].iloc[0],
    bf_params['BETA'].iloc[0],
    bf_params['Ks'].iloc[0],
    bf_params['Kb'].iloc[0],
    bf_params['Kz'].iloc[0]
]

flow = [
    bf_params['Qthresh'].iloc[0],
    bf_params['Rs'].iloc[0],
    bf_params['Rb1'].iloc[0],
    bf_params['Rb2'].iloc[0],
    bf_params['Prec'].iloc[0],
    bf_params['Frac4Rise'].iloc[0]
]

# Use parameters for baseflow separation
SBT = bfs.base_table(basin_char[1], basin_char[2], basin_char[3],
                     gw_hyd[1], gw_hyd[3], streamflow_data, basin_char[4])
result = bfs.PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
```

## Tips and Best Practices

- **Data Requirements**: Ensure your streamflow data has columns named 'Date' and 'Streamflow'
- **Parameter Files**: Site parameter files should contain basin characteristics, groundwater hydraulic parameters, and flow metrics
- **Calibration Period**: Use a representative period for calibration that captures different flow conditions (at least several months, preferably a full year)
- **Calibration Data Quality**: Ensure high-quality streamflow data with minimal gaps for best calibration results
- **Calibration Time**: Calibration can take several minutes to hours depending on data length; be patient
- **Initial Conditions**: For forecasting, always extract initial conditions from the last time step of your calibration run
- **Parameter Validation**: After calibration, check that parameters are physically reasonable (e.g., conductivities between 10⁻⁸ and 10⁵ m/day)
- **Error Metrics**: Review the calibration error and BFF values to assess calibration quality
- **Visualization**: Use the provided plotting functions to visualize and validate your results
- **Saving Parameters**: Save calibrated parameters for reuse to avoid re-calibrating for the same site