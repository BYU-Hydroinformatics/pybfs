# PyBFS Documentation

Welcome to PyBFS - a Python implementation of the USGS Baseflow Separation (BFS) model for hydrological analysis and forecasting.

## Overview

PyBFS is a Python implementation of the BFS (Baseflow Separation) model, a non-linear state-space model originally developed by the U.S. Geological Survey (USGS). The model provides a physically-based approach to separating streamflow into its component parts and forecasting baseflow during dry periods.

### Understanding Baseflow Separation

Baseflow represents the steady, sustained component of streamflow that is primarily contributed by groundwater discharge to streams and rivers. Understanding and quantifying baseflow is essential for managing water resources during drought conditions, assessing groundwater-surface water interactions, estimating groundwater residence times, forecasting low-flow conditions, and supporting ecological flow assessments. Traditional approaches to baseflow separation have relied on graphical methods that are often subjective and inadequate for complex hydrological environments such as snowmelt-dominated rivers or watersheds with complicated geology. The BFS model addresses these limitations through a sophisticated numerical approach grounded in physical principles.

![figure_1.png](images/figure_1.png)

### Baseflow Forecasting with PyBFS

One of the most powerful features of PyBFS is its ability to forecast future baseflow without requiring precipitation data. This capability is particularly valuable during drought conditions or for planning purposes when you need to understand how streamflow will evolve in the absence of significant rainfall events. The forecasting workflow begins by running the baseflow separation algorithm on a historical calibration period to establish the current state of the watershed's surface and subsurface reservoirs. From the final time step of this calibration period, the model extracts a complete snapshot of the system state including water levels in both reservoirs, storage volumes, flow rates, and recharge conditions.

Using these initial conditions, the forecast function projects the system forward in time by simulating how the reservoirs will drain and interact according to the established physical relationships. The model assumes zero precipitation during the forecast period, making it essentially a recession forecast that shows how baseflow will decline based purely on the drainage of existing storage. This approach is particularly useful for answering questions like "How long will this stream maintain a minimum flow threshold?" or "What will baseflow conditions be like in 30, 60, or 90 days if no significant rain occurs?" Water resource managers can use these forecasts to make informed decisions about water allocation, flow augmentation needs, or drought emergency declarations.

![figure_3.png](images/figure_3.png)

### The BFS Model Approach

The BFS model employs a **non-linear state-space framework** in which baseflow is conceptualized as a function of upstream storage—an unmeasured state variable that the model infers from observed streamflow patterns. This approach represents a significant advancement over traditional methods because it captures the complex, non-linear relationship between groundwater storage and baseflow discharge. The model uses a coupled reservoir system where surface and subsurface reservoirs interact through infiltration and recharge processes, creating a physically realistic representation of watershed behavior. Model parameters directly relate to measurable basin characteristics such as geometry and hydraulic properties, providing a physical basis for the separation rather than relying on arbitrary assumptions.

The optimization framework maximizes the baseflow component while maintaining physical constraints, ensuring that the separated components are both mathematically optimal and physically plausible. A particularly valuable feature is the model's ability to forecast baseflow into future periods without requiring precipitation data, making it especially useful for drought analysis and low-flow forecasting. The model partitions total observed streamflow into three distinct components: **baseflow** (the groundwater contribution that is sustained and slowly varying), **surface flow** (water moving as overland flow and through shallow subsurface pathways), and **direct runoff** (the quick response to precipitation events).

### Original USGS Research

PyBFS implements the methodology described in the USGS Scientific Investigations Report by Christopher Konrad: **BFS—A non-linear state-space model for baseflow separation and prediction** (Report 2022–5114, published December 2022, [https://doi.org/10.3133/sir20225114](https://doi.org/10.3133/sir20225114)). The original BFS model was developed at the USGS and implemented in the R programming language. What makes this model particularly robust is that it was calibrated and tested across an extensive network of 13,208 USGS streamgages distributed throughout the United States, spanning diverse hydrological contexts from arid regions to humid climates, and from groundwater-dominated systems to snowmelt-driven rivers. This extensive validation provides confidence that the model can perform reliably across a wide range of environmental conditions.

### PyBFS Implementation

This Python implementation brings the BFS methodology to the Python scientific computing ecosystem. The package includes the core BFS functions for baseflow table generation, the separation algorithm itself, and forecasting capabilities. Built-in visualization tools allow users to quickly plot and analyze their results without needing to write custom plotting code. The implementation includes utilities for parameter management that help users extract and organize the site-specific parameters required by the model. To enable efficient computation even with long time series, the implementation uses a lookup table approach that allows rapid baseflow calculations without repeatedly solving complex equations at each time step.

## Quick Links

[Installation](user-guide/installation.md) | [User Guide](user-guide/overview.md) | [Usage Examples](user-guide/examples.md) | [API Reference](api/reference.md)

## Key Features

PyBFS is a **physically-based model** that uses actual basin geometry and hydraulic properties rather than empirical coefficients. The **non-linear state-space framework** captures the complex dynamics of baseflow response to changing storage conditions. The **baseflow forecasting** capability allows projection of future baseflow during dry periods when precipitation is minimal or absent. The model provides **comprehensive output** including not just baseflow but also surface flow, direct runoff, and internal storage states, giving users a complete picture of watershed hydrology. Built-in **visualization functions** simplify the process of examining model results and comparing simulated versus observed streamflow. The approach is **well-tested**, having been calibrated and validated across thousands of USGS streamgages representing diverse hydrological settings.

## Applications

PyBFS is particularly valuable for hydrological research and watershed analysis where understanding the sources of streamflow is critical. Water resources managers use the model for planning and management activities, especially during drought conditions when baseflow becomes the dominant or only source of streamflow. The forecasting capability makes it useful for drought monitoring and low-flow forecasting applications. Researchers studying groundwater-surface water interactions can use the model to quantify the exchange between these systems. Environmental flow assessments benefit from the model's ability to separate the sustained baseflow component that is often critical for aquatic ecosystems. Climate change impact studies can use the model to understand how changing conditions might affect baseflow contributions over time.

## Getting Help

For questions, issues, or contributions, please visit the [GitHub repository](https://github.com/BYU-Hydroinformatics/pybfs).

## Citation

If you use PyBFS in your research, please cite both this implementation and the original USGS publication:

> Konrad, C.P., 2022, BFS—A non-linear state-space model for baseflow separation and prediction: U.S. Geological Survey Scientific Investigations Report 2022–5114, 28 p., https://doi.org/10.3133/sir20225114