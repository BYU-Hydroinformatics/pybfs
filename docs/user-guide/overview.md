# User Guide Overview

This section provides comprehensive documentation on how to use PyBFS for baseflow separation and streamflow forecasting.

## Introduction

PyBFS (Python Baseflow Separation) is a physically-based hydrological model for separating streamflow into its component parts and forecasting future baseflow. The model uses a coupled surface-subsurface reservoir approach to partition total streamflow into:

- **Baseflow**: Water contributed by groundwater discharge
- **Surface Flow**: Water flowing across the land surface
- **Direct Runoff**: Quick response flow from precipitation events

## Core Concepts

### Two-Reservoir System

PyBFS uses a dual-reservoir conceptual model:

1. **Surface Reservoir**: Represents water stored on and near the land surface, including soil moisture and shallow subsurface storage
2. **Base Reservoir**: Represents deeper groundwater storage that contributes to baseflow

These reservoirs interact through infiltration (surface to base) and recharge processes, creating a physically-realistic representation of watershed hydrology.

### Baseflow Table

The baseflow table (SBT) is a lookup table that relates groundwater levels to storage and discharge. It is generated using basin geometry and hydraulic properties:

- Basin length (lb) and width (wb)
- Aquifer shape parameter (beta)
- Base hydraulic conductivity (kb)
- Porosity (por)

The table enables efficient computation of baseflow for different groundwater conditions without solving complex groundwater equations at each time step.

### Parameter Sets

PyBFS requires three parameter sets:

**1. Basin Characteristics** (`basin_char`):

- Area: Basin drainage area (m²)
- lb: Basin length (m)
- x1: Initial longitudinal position (m)
- wb: Basin width (m)
- por: Porosity (0-1)

**2. Groundwater Hydraulic Parameters** (`gw_hyd`):

- alpha: Surface reservoir shape parameter
- beta: Base reservoir shape parameter
- ks: Surface hydraulic conductivity (m/day)
- kb: Base hydraulic conductivity (m/day)
- kz: Vertical hydraulic conductivity (m/day)

**3. Flow Metrics** (`flow`):

- qthresh: Flow threshold for recession identification
- rs: Recession slope parameter
- rb1, rb2: Baseflow recession parameters
- prec: Precision threshold for convergence
- fr4rise: Fraction threshold for rise detection

### Model Workflow

The typical PyBFS workflow consists of:

1. **Load Data**: Import streamflow observations and site parameters
2. **Extract Parameters**: Use `get_values_for_site()` to extract parameters for your site
3. **Generate Baseflow Table**: Create the lookup table with `base_table()`
4. **Run Separation**: Execute `PyBFS()` to separate streamflow components
5. **Visualize**: Plot results using provided plotting functions
6. **Forecast** (optional): Use `forecast()` to predict future baseflow

### Forecasting

PyBFS can forecast future baseflow given initial conditions. The forecasting approach:

- Requires a calibration period to establish initial reservoir states
- Extracts initial conditions (water levels, storage, flows) from the last calibration time step
- Projects baseflow forward in time based on the reservoir dynamics
- Does not require future precipitation data (uses zero precipitation assumption)

This makes it useful for drought analysis and low-flow forecasting.

## Best Practices

- **Data Quality**: Ensure your streamflow data is complete and quality-checked before analysis
- **Parameter Estimation**: Site parameters should be based on measured basin characteristics when possible
- **Calibration Length**: Use a calibration period of at least several months to capture seasonal variations
- **Validation**: Always validate model results against observed data for a holdout period
- **Physical Realism**: Check that separated components (baseflow, surface flow) are physically reasonable
- **Time Series Continuity**: Ensure your streamflow data has no gaps during the analysis period
- **Units Consistency**: Maintain consistent units throughout (meters and days are standard)