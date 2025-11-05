# Usage Examples

This page demonstrates practical applications of PyBFS for baseflow separation and forecasting.

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

## Google Colab Example

You can also run PyBFS in a Google Colab environment. Here is an example notebook that demonstrates how to do this:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/BYU-Hydroinformatics/pybfs/blob/main/notebooks/pybfs_sample.ipynb){target="_blank"}

Before running this notebook, make sure to upload the necessary data files to the Colab environment. First of all, you will need to upload a copy of the pybfs.py file to the Colab environment. You can do this by clicking on the folder icon on the left sidebar, then clicking the upload icon (a paper with an upward arrow) to upload the pybfs.py file from your local machine.

The pybfs.py file can be found in the main PyBFS GitHub repository here: [pybfs.py](https://github.com/BYU-Hydroinformatics/pybfs/blob/main/pybfs.py)

Next, you will need to upload the data files used in the examples above (e.g., `2312200_data.csv`, `bfs_params_50.csv`). You can upload these files in the same way you uploaded the pybfs.py file. You can download a copy of the files using these links:

[2312200_data.csv](../files/2312200_data.csv)<br>
[bfs_params_50.csv](../files/bfs_params_50.csv)

Once you have uploaded the necessary files, you can run the cells in the notebook to perform baseflow separation and forecasting using PyBFS.

## Tips and Best Practices

- **Data Requirements**: Ensure your streamflow data has columns named 'Date' and 'Streamflow'
- **Parameter Files**: Site parameter files should contain basin characteristics, groundwater hydraulic parameters, and flow metrics
- **Calibration Period**: Use a representative period for calibration that captures different flow conditions
- **Initial Conditions**: For forecasting, always extract initial conditions from the last time step of your calibration run
- **Visualization**: Use the provided plotting functions to visualize and validate your results