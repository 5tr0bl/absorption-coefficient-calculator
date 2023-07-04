import numpy as np
from numpy import cos, sin, pi, sqrt, abs, dot
import pandas as pd
import plotly.graph_objects as go
import pendulum
import streamlit as st

# Define a class that works with three variables and returns a dataframe and a plot
class absorption_coefficient:
    """This class takes three variables as input and returns a plotly plot and a dataframe.
    """

    def __init__(
            self,
            f_start:int, 
            f_end:int, 
            d:int):
        self.f_start = f_start
        self.f_end = f_end
        self.d = d

    def plot(self):
        fig = go.Figure()
        x = np.arange(self.f_start, self.f_end + 10, 10)
        alpha = 0.9 * np.log10(x*self.d) - 2.4
        alpha[alpha > 1] = 1
        alpha[alpha < 0] = 0
        global df
        df = pd.DataFrame({'Frequency': x, 'Absorption Coefficient': alpha})
        fig.add_trace(go.Scatter(x=x, y=alpha, mode='lines'))
        fig.update_layout(yaxis_range=[0, 1.1])
        fig.update_xaxes(showgrid=True)
        fig.update_layout(title='Absorption coefficient of a {} mm material'.format(self.d),
                          xaxis_title='Frequency in [Hz]',
                          yaxis_title='Absorption Coefficient')
        return fig
    
    def data(self):
        return df


# Define a function that converts a dataframe to a CSV file
@st.cache(show_spinner=False)
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

# Define a function that creates an input field for a number
def create_input_field():
    """Creates a Streamlit input field to enter a number.

    Returns:
        int: Streamlit input field returning an integer.
    """
    return st.number_input(
        label="Select material Nr. 2 thickness [mm]:",
        min_value=0,
        max_value=100,
        value=50,
        step=1,
    )

# Define the transfer matrix method function
def tmm(
        theta: float, # needs to be in radians for internal calculation
        k0: float,
        k1: float,
        L1: float,
        Z1: float,
        k2: float,
        L2: float,
        Z2: float
    ) -> np.array:
    
    ## neue k_pz bei Eingallswinkel theta
    k_0x = k0*sin(theta)
    k_1x = k_0x
    # k_2x = k2*sin(theta)

    k_1z = sqrt(k1**2 - k_1x**2)
    # k_2z = sqrt(k2**2 - k_2x**2)

    # T1, T2, T3 definieren

    T1 = np.array([[cos(k_1z*L1), 1j*Z1*(k1/k_1z)*sin(k_1z*L1)],
                   [(1j/Z1)*(k_1z/k1)*sin(k_1z*L1), cos(k_1z*L1)]])

    '''
    T2 = np.array([[cos(k_2z*L2), 1j*Z2*(k2/k_2z)*sin(k_2z*L2)],
                   [(1j/Z2)*(k_2z/k2)*sin(k_2z*L2), cos(k_2z*L2)]])
    '''


    # Rigid Backed
    T3 = np.array([[1, 0],  
                   [0, 1]])

    # T_total berechnen
    T_total = np.dot(T1, T2)

    return T_total

def plotly_go_line(x,y,x_label,y_label,title):
    """Creates a plotly-go line plot.

    Returns:
        plotly.graph_objects.Figure: Plotly line plot.
    """
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines'))
    fig.update_xaxes(showgrid=True)
    fig.update_layout(title=title,
                    xaxis_title=x_label,
                    yaxis_title=y_label,
                    width=1000,
                    height=500)
    return fig