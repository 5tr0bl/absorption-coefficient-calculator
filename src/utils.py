import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pendulum
import streamlit as st


# Define a function that converts a dataframe to a CSV file
# @st.cache_data(show_spinner=False)
def _convert_df(df: pd.DataFrame):
    return df.to_csv().encode("utf-8")


# Define a function that creates a download button for a dataframe
def create_df_export_button(
        df: pd.DataFrame,
        title: str,
        ts: pendulum.DateTime | None,
) -> bool:
    """Creates a Streamlit button to export a dataframe to a CSV file.

    Args:
        df (pd.DataFrame): Dataframe to export.
        title (str): Title of the data.
        ts (pendulum.DateTime): Optional datetime that will be used for file name.

    Returns:
        bool: Streamlit button functioning a boolean type.
    """

    if ts is None:
        ts = pendulum.now()

    ts_formatted = ts.to_datetime_string().translate(
        str.maketrans(
            {
                "/": "-",
                ":": "-",
            }
        )
    )

    file_name = f"{title}_{ts_formatted}.csv".replace(" ", "_").lower()

    return st.download_button(
        label="Export",
        data=_convert_df(df=df),
        file_name=file_name,
        mime="text/csv",
    )


# Define transfer matrix methode
def tmm(
        theta: float,
        k0: float,
        k1: float,
        L1: float,
        Z1: float,
        k2: float,
        L2: float,
        Z2: float
):
    ## neue k_pz bei Eingallswinkel theta
    k_0x = k0 * np.sin(theta)
    k_1x = k_0x
    k_2x = k2 * np.sin(theta)

    k_1z = np.sqrt(k1 ** 2 - k_1x ** 2)
    k_2z = np.sqrt(k2 ** 2 - k_2x ** 2)

    # T1, T2, T3 definieren

    T1 = np.array([[np.cos(k_1z * L1), 1j * Z1 * (k1 / k_1z) * np.sin(k_1z * L1)],
                   [(1j / Z1) * (k_1z / k1) * np.sin(k_1z * L1), np.cos(k_1z * L1)]])

    T2 = np.array([[np.cos(k_2z * L2), 1j * Z2 * (k2 / k_2z) * np.sin(k_2z * L2)],
                   [(1j / Z2) * (k_2z / k2) * np.sin(k_2z * L2), np.cos(k_2z * L2)]])

    # Rigid Backed
    T3 = np.array([[1, 0],
                   [0, 1]])

    # T_total berechnen
    T_total = np.dot(T1, T2)

    return T_total


def plotly_go_line(x, y, x_label, y_label, title):
    """Creates a plotly-go line plot.

    Returns:
        plotly.graph_objects.Figure: Plotly line plot.
    """
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines'))
    fig.update_xaxes(showgrid=True, type="log")
    fig.update_layout(title=title,
                      xaxis_title=x_label,
                      yaxis_title=y_label,
                      yaxis_range=[0, 1],
                      width=1000,
                      height=500)
    return fig


def plotly_freq_bands(x, y, x_label, y_label, title, plot_type='oct'):
    """Creates a plotly-go bar plot for octave bands.

    Returns:
        plotly.graph_objects.Figure: Plotly bar plot.
    """
    # mapped_alphas = dict(zip(x,y))

    if plot_type == 'oct':
        center_freqs = [31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]
        bw_factor = 2
    elif plot_type == 'third':
        center_freqs = [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200,
                        250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
                        2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000]
        bw_factor = 3
    else:
        raise ValueError("Invalid Plot Type")

    freq_bands = []

    for center_freq in center_freqs:
        lower_cutoff = center_freq / bw_factor
        upper_cutoff = center_freq * bw_factor

        freq_band = {
            "center_frequency": center_freq,
            "lower_cutoff_frequency": lower_cutoff,
            "upper_cutoff_frequency": upper_cutoff
        }

        freq_bands.append(freq_band)

    alphas_mean = []  # will hold the mean value for all alphas in each freq band

    for band in freq_bands:
        band_y_values = y[(x >= band['lower_cutoff_frequency']) & (x <= band['upper_cutoff_frequency'])]
        mean_y_value = np.mean(band_y_values) if len(band_y_values) > 0 else 0
        alphas_mean.append(mean_y_value)

    # Create evenly spaced x-axis values for plotting
    x_ticks = np.arange(len(center_freqs))

    # Create bar plot
    fig = go.Figure(data=go.Bar(x=x_ticks, y=alphas_mean))

    # Set the x-axis tick positions and labels
    fig.update_layout(
        xaxis=dict(
            tickvals=x_ticks,
            ticktext=center_freqs,
            title=x_label
        ),
        yaxis=dict(
            range=[0, 1],
            title=y_label
        ),
        title=title
    )

    return fig
