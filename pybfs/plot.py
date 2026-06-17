# -*- coding: utf-8 -*-
"""Plotting functions for PyBFS

Functions for visualizing baseflow separation results and forecasts.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


def plot_baseflow_simulation(
    streamflow,
    baseflow,
    title="Baseflow Simulation",
    figsize=(10, 6),
    save=True,
    dpi=300,
):
    """Plots observed streamflow vs simulated baseflow from PyBFS

    Creates a time series plot comparing observed total streamflow with simulated
    baseflow. Useful for visualizing baseflow separation results and assessing
    model performance.

    Parameters
    ----------
    streamflow : pd.DataFrame
        DataFrame containing observed data with columns Date (datetime) and Streamflow
        (observed streamflow m³/day)
    baseflow : pd.DataFrame
        Output from PyBFS() containing flow components with column Baseflow (simulated
        baseflow component m³/day)
    title : str, optional
        Plot title, default is "Baseflow Simulation"


    Notes
    -----
    Streamflow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Black line shows total observed streamflow. Green line shows simulated baseflow component.
    Displays plot using matplotlib.

    Examples
    --------
    >>> results = bfs(streamflow_data, SBT, basin_char, gw_hyd, flow)
    >>> plot_df = plot_baseflow_simulation(streamflow_data, results)
    """
    # Prepare DataFrame for plotting
    # df = pd.DataFrame({
    #     "date": pd.to_datetime(streamflow["Date"]),
    #     "streamflow": streamflow["Streamflow"],
    #     "baseflow": tmp["Baseflow"]
    # })

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    # Plot observed streamflow (converted from m³/day to m³/s)
    ax.plot(streamflow["Date"], streamflow["Streamflow"] / 86400, color="black",
            label="Streamflow", linewidth=1)

    # Plot simulated baseflow (converted similarly)
    ax.plot(streamflow["Date"], baseflow["Baseflow"] / 86400, color="green",
            label="PyBFS", linewidth=1.5)

    # Labels and formatting
    ax.set_xlabel("Date")
    ax.set_ylabel("Flow (cms)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.tick_params(axis="both", which="major")

    # Date formatting
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    fig.autofmt_xdate()

    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    if save:
        fig.savefig("baseflow_simulation.png", dpi=dpi, bbox_inches="tight")
    plt.show()

def plot_forecast(
    training_streamflow,
    training_bfs,
    forecast_data,
    title="Baseflow Training + Forecast",
    streamflow_col="Streamflow",
    figsize=(10, 6),
    save=True,
    dpi=300,
):
    """Plot training-period streamflow/baseflow and forecast-period baseflow.

    This helper is designed to match the workflow in `main.py`:
    - Run `bfs()` on a training/calibration period to obtain baseflow separation
    - Run `forecast()` on a future period using initial conditions from the training run
    - Plot training streamflow + training baseflow (solid lines), and forecast baseflow
      for the forecast period (dashed line) on the same axis.

    Parameters
    ----------
    training_streamflow : pd.DataFrame
        Observed streamflow for the training period with columns:
        - 'Date': datetime-like
        - `streamflow_col`: observed streamflow (m³/day)
    training_bfs : pd.DataFrame
        Output from `bfs()` for the training period with columns:
        - 'Date': datetime-like
        - 'Baseflow': simulated baseflow (m³/day)
    forecast_data : pd.DataFrame
        Output from `forecast()` with columns:
        - 'Date': datetime-like
        - 'Baseflow': forecast baseflow (m³/day)
    title : str, optional
        Plot title.
    streamflow_col : str, optional
        Column name for observed streamflow in `training_streamflow`.


    Notes
    -----
    Values are converted from m³/day to cms (m³/s) by dividing by 86400.
    """
    # Training period merge (align on Date)
    train_stream = pd.DataFrame(
        {
            "Date": pd.to_datetime(training_streamflow["Date"]),
            "Streamflow": training_streamflow[streamflow_col],
        }
    )
    train_bf = pd.DataFrame(
        {
            "Date": pd.to_datetime(training_bfs["Date"]),
            "Baseflow": training_bfs["Baseflow"],
        }
    )
    train_df = (
        pd.merge(train_stream, train_bf, on="Date", how="inner")
        .sort_values("Date")
        .reset_index(drop=True)
    )

    # Forecast period
    f_df = pd.DataFrame(
        {
            "Date": pd.to_datetime(forecast_data["Date"]),
            "Baseflow": forecast_data["Baseflow"],
        }
    ).sort_values("Date")

    fig, ax = plt.subplots(figsize=figsize)

    # Training: observed streamflow + simulated baseflow (solid)
    ax.plot(
        train_df["Date"],
        train_df["Streamflow"] / 86400,
        color="black",
        label="Streamflow (training)",
        linewidth=1,
    )
    ax.plot(
        train_df["Date"],
        train_df["Baseflow"] / 86400,
        color="green",
        label="Baseflow (training)",
        linewidth=1.5,
    )

    # Forecast: baseflow (dashed)
    ax.plot(
        f_df["Date"],
        f_df["Baseflow"] / 86400,
        color="green",
        linestyle="--",
        label="Baseflow (forecast)",
        linewidth=1.5,
    )

    # Visual separator at forecast start
    if len(f_df) > 0:
        ax.axvline(f_df["Date"].iloc[0], color="gray", linestyle=":", linewidth=1)

    ax.set_xlabel("Date")
    ax.set_ylabel("Flow (cms)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.tick_params(axis="both", which="major")

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    fig.autofmt_xdate()

    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    if save:
        fig.savefig("figure_3.png", dpi=dpi, bbox_inches="tight")
    plt.show()


def plot_forecast_baseflow(forecast_data, figsize=(9, 5), save=True, dpi=300):
    """Plot baseflow forecast time series

    Creates a time series plot of forecasted baseflow. Shows how baseflow
    is expected to evolve over the forecast period based on drainage from
    storage reservoirs.

    Parameters
    ----------
    forecast_data : pd.DataFrame
        Forecast results from forecast() function with columns Date (datetime of forecast)
        and Baseflow (forecasted baseflow m³/day)

    Notes
    -----
    Baseflow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Green line shows forecasted baseflow. Displays plot using matplotlib. No observed
    streamflow is shown (forecast period has no observations).

    Examples
    --------
    >>> forecast_result = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
    >>> plot_forecast_baseflow(forecast_result)
    """
    fig, axs = plt.subplots(figsize=figsize)
    date = pd.to_datetime(forecast_data["Date"])
    axs.plot(date, forecast_data['Baseflow']/86400, color='green', label='PyBFS Baseflow', linewidth=1.5)

    # Add legend
    axs.legend(loc='upper right')

    # Set title and axis labels
    axs.set_title(f"Baseflow Forecast")
    axs.set_xlabel('Date')
    axs.set_ylabel('Flow (cms)')

    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()

    # Tick label font sizes
    axs.tick_params(axis='x')
    axs.tick_params(axis='y')

    if save:
        fig.savefig("baseflow_forecast.png", dpi=dpi, bbox_inches="tight")
    plt.show()


def plot_forecast_baseflow_streamflow(
    forecast_data, streamflow, figsize=(9, 5), save=True, dpi=300
):
    """Plot forecast baseflow with observed streamflow for comparison

    Creates a time series plot comparing forecasted baseflow against observed
    streamflow. Useful for validating forecast performance when observations
    are available for the forecast period.

    Parameters
    ----------
    forecast_data : pd.DataFrame
        Forecast results from forecast() function with columns Date (datetime of forecast)
        and Baseflow (forecasted baseflow m³/day)
    streamflow : pd.DataFrame
        Observed streamflow data for the forecast period with columns Date (datetime of
        observations) and Streamflow (observed streamflow m³/day)

    Notes
    -----
    Flow values are converted from m³/day to m³/s by dividing by 86400 seconds/day.
    Blue line shows observed total streamflow. Green line shows forecasted baseflow.
    Displays plot using matplotlib. Used for forecast validation when observations become available.

    Examples
    --------
    >>> forecast_result = forecast(forecast_df, SBT, basin_char, gw_hyd, flow, ini)
    >>> # Get observed data for same period
    >>> obs_data = streamflow_data[(streamflow_data['Date'] >= '2018-10-01') &
    ...                             (streamflow_data['Date'] <= '2018-11-30')]
    >>> plot_forecast_baseflow_streamflow(forecast_result, obs_data)
    """
    fig, axs = plt.subplots(figsize=figsize)
    date = pd.to_datetime(forecast_data["Date"])

    axs.plot(date, streamflow['Streamflow']/86400, color='blue', label='USGS Streamflow', linewidth=1.5)
    axs.plot(date, forecast_data['Baseflow']/86400, color='green', label='PyBFS Baseflow', linewidth=1.5)

    # Add legend
    axs.legend(loc='upper right')

    # Set title and axis labels
    axs.set_title(f"Baseflow Forecast")
    axs.set_xlabel('Date')
    axs.set_ylabel('Flow (cms)')

    axs.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()

    # Tick label font sizes
    axs.tick_params(axis='x')
    axs.tick_params(axis='y')
    #end

    if save:
        fig.savefig("baseflow_forecast_streamflow.png", dpi=dpi, bbox_inches="tight")
    plt.show()

