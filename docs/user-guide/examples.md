# Usage Examples

This page demonstrates practical applications of PyBFS for baseflow separation and forecasting.

## Complete Workflow Example

This example demonstrates the full PyBFS workflow including:

1. Loading streamflow data and site parameters
2. Generating a baseflow table
3. Running PyBFS for baseflow separation
4. Visualizing results
5. Creating forecasts

```python
#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pybfs import (
    get_values_for_site,
    base_table,
    PyBFS,
    plot_baseflow_simulation,
    forecast,
    plot_forecast_baseflow,
    plot_forecast_baseflow_streamflow
)

# Load streamflow data
streamflow_data = pd.read_csv('tutorial/2312200_data.csv')

# Load site parameters
bfs_params_usgs = pd.read_csv('tutorial/bfs_params_50.csv')

# Get parameters for specific site
site_number = 2312200
basin_char, gw_hyd, flow = get_values_for_site(bfs_params_usgs, site_number)

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
SBT = base_table(lb, x1, wb, beta, kb, streamflow_data, por)
print(f"Baseflow table generated with {len(SBT)} rows")

# Run PyBFS
result = PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)
print(f"PyBFS completed for {len(result)} time steps")

# Display summary statistics
print("\n=== Results Summary ===")
print(f"Total observed flow: {result['Qob'].sum():.2f}")
print(f"Total simulated flow: {result['Qsim'].sum():.2f}")
print(f"Total baseflow: {result['Baseflow'].sum():.2f}")
print(f"Total surface flow: {result['SurfaceFlow'].sum():.2f}")
print(f"Total direct runoff: {result['DirectRunoff'].sum():.2f}")

# Plot results
plot_baseflow_simulation(streamflow_data, result)
```

## Forecasting Example

Once you have calibrated the model, you can create forecasts for future periods:

```python
# Filter data for calibration period (Jan-Sep 2018)
start_date = '2018-01-01'
end_date = '2018-09-30'
streamflow_data['Date'] = pd.to_datetime(streamflow_data['Date'])
streamflow_data_filtered = streamflow_data[
    (streamflow_data['Date'] >= start_date) & (streamflow_data['Date'] <= end_date)
]

# Run PyBFS for calibration period
tmp2 = PyBFS(streamflow_data_filtered, SBT, basin_char, gw_hyd, flow)

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
f = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
print(f"Forecast completed for {len(f)} time steps")

# Plot forecast
plot_forecast_baseflow(f)

# Plot forecast with observed data for comparison
forecast_start = '2018-10-01'
forecast_end = '2018-11-30'
streamflow_data_forecast = streamflow_data[
    (streamflow_data['Date'] >= forecast_start) & (streamflow_data['Date'] <= forecast_end)
]

plot_forecast_baseflow_streamflow(f, streamflow_data_forecast)
```

## Tips and Best Practices

- **Data Requirements**: Ensure your streamflow data has columns named 'Date' and 'Streamflow'
- **Parameter Files**: Site parameter files should contain basin characteristics, groundwater hydraulic parameters, and flow metrics
- **Calibration Period**: Use a representative period for calibration that captures different flow conditions
- **Initial Conditions**: For forecasting, always extract initial conditions from the last time step of your calibration run
- **Visualization**: Use the provided plotting functions to visualize and validate your results