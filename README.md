# PyBFS

PyBFS - Python Baseflow Separation

A Python implementation of physically-based Baseflow Separation algorithms for hydrological analysis. PyBFS uses a coupled surface-subsurface reservoir model to separate total streamflow into surface flow, baseflow, and direct runoff components.

## Features

- Physically-based baseflow separation using coupled reservoir models
- Forecast baseflow using initial conditions
- Comprehensive visualization tools
- Well-documented API with NumPy-style docstrings

## Installation

### From PyPI (recommended)

```bash
pip install pybfs
```

### From source

```bash
git clone https://github.com/BYU-Hydroinformatics/pybfs.git
cd pybfs
pip install .
```

### Development installation

```bash
git clone https://github.com/BYU-Hydroinformatics/pybfs.git
cd pybfs
pip install -e .
```

## Quick Start

```python
import pandas as pd
from pybfs import get_values_for_site, base_table, PyBFS, plot_baseflow_simulation

# Load streamflow data
streamflow_data = pd.read_csv('streamflow.csv')

# Load site parameters
params_df = pd.read_csv('site_parameters.csv')

# Get parameters for specific site
site_number = 2312200
basin_char, gw_hyd, flow = get_values_for_site(params_df, site_number)

# Extract parameters
area, lb, x1, wb, por = basin_char[0], basin_char[1], basin_char[2], basin_char[3], basin_char[4]
alpha, beta, ks, kb, kz = gw_hyd[0], gw_hyd[1], gw_hyd[2], gw_hyd[3], gw_hyd[4]

# Generate baseflow table
SBT = base_table(lb, x1, wb, beta, kb, streamflow_data, por)

# Run PyBFS baseflow separation
result = PyBFS(streamflow_data, SBT, basin_char, gw_hyd, flow)

# Visualize results
plot_baseflow_simulation(streamflow_data, result)
```

## Requirements

- Python 3.7 or higher
- numpy >= 1.19.0
- pandas >= 1.1.0
- matplotlib >= 3.3.0

## Documentation

Full documentation is available at: [Documentation](https://github.com/BYU-Hydroinformatics/pybfs)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/contributing.md) for guidelines.

## Citation

If you use PyBFS in your research, please cite:

```
PyBFS - Python Baseflow Separation
BYU Hydroinformatics
https://github.com/BYU-Hydroinformatics/pybfs
```

## Support

For issues, questions, or contributions, please visit: https://github.com/BYU-Hydroinformatics/pybfs/issues
