# -*- coding: utf-8 -*-
"""PyBFS - Python Baseflow Separation

A Python implementation of Baseflow Separation algorithms for hydrological analysis.
"""

__version__ = "0.1.1"
__author__ = "BYU Hydroinformatics"
__email__ = ""

# Import from utilities module
from .utilities import (
    sur_z,
    sur_store,
    sur_q,
    dir_q,
    infiltration,
    recharge,
    get_values_for_site,
    base_table,
    forecast,
    flow_metrics,
    bf_ci,
)

# Import from bfs module
from .bfs import bfs

# Import from plot module
from .plot import (
    plot_baseflow_simulation,
    plot_forecast,
    plot_forecast_baseflow,
    plot_forecast_baseflow_streamflow,
)

# Import from skill module
from .skill import (
    modified_strict_baseflow,
    separation_skill,
    forecast_skill,
)

# Import from calibrate module
from .calibrate import (
    bfs_calibrate,
    objective,
    ini_params,
    cal_initial,
    cal_basetable,
    cal_base,
    cal_surface,
)

__all__ = [
    "sur_z",
    "sur_store",
    "sur_q",
    "dir_q",
    "infiltration",
    "recharge",
    "get_values_for_site",
    "base_table",
    "flow_metrics",
    "bf_ci",
    "bfs",
    "forecast",
    "bfs_calibrate",
    "objective",
    "ini_params",
    "cal_initial",
    "cal_basetable",
    "cal_base",
    "cal_surface",
    "modified_strict_baseflow",
    "separation_skill",
    "forecast_skill",
    "plot_baseflow_simulation",
    "plot_forecast",
    "plot_forecast_baseflow",
    "plot_forecast_baseflow_streamflow",
]
