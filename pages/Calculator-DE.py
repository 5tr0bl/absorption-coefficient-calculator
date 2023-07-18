import streamlit as st
import numpy as np
import pandas as pd

from src import utils, models, absorptioncoeff

st.set_page_config(
    page_title="Absorptionsgrad Rechner DE",
    layout="wide"
)
# --- Initialising SessionState ---
if "load_state" not in st.session_state:
    st.session_state.load_state = False

# Set the title and logo of the app
col1, col2 = st.columns(2)
col1.title('Absorptionsgrad Rechner')
col2.image('img/logo-tuberlin-header2.png', width=350)

################## Input Section ##################
# st.markdown('----')

# Define the dropdown menus for the frequency
st.header('Globale Parameter :globe_with_meridians:')
with st.expander('Werte ein/ausblenden...'):
    st.markdown('##### Frequenzbereich')
    col1, col2 = st.columns(2)
    f_min, f_max = col1.slider('Anfangs- und Endfrequenz [Hz]', 0, 20000, (0, 10000), step=10)
    f_range = np.arange(f_min, f_max, 1)
    f_range_full = np.arange(1, 20000, 1)
    plot_type = col2.selectbox('Darstellung', ('Graph', 'Oktavbänder', 'Terzbänder'))

    col1, col2, col3 = st.columns(3)
    col1.markdown('##### Lufttemperatur')
    air_temp = col1.number_input('in [°C]', step=1, value=20)
    col2.markdown('##### Luftdruck')
    air_pressure = col2.number_input('in [Pa]', step=1, value=101325)
    col3.markdown('##### Einfallswinkel')
    theta = col3.number_input('in [°]', step=1, value=0)

st.markdown('----')

# Materialen Eingabe
st.header('Material Parameter :hammer:')
with st.expander('Werte ein/ausblenden...'):
    col1, col2, col3 = st.columns(3)
    num_materials = col1.number_input("Wähl die Anzahl an Materialen:", min_value=1, max_value=5, value=2, step=1)
    material_dict = {}

    # Create the specified number of columns for each material
    columns = st.columns(num_materials)
    # Display variables in each column and save them
    for i, column in enumerate(columns):
        column.markdown(f"##### Material # {i + 1}")
        key = f"Material {i + 1}"
        value1 = column.selectbox(f"Modell", options=['Bitte wählen Sie', 'Poröser', 'Lochplatte', 'Platte', 'Luft'],
                                  key=f"value1_{i}")
        value2 = column.number_input(f"Dicke [mm]", key=f"value2_{i}", format='%0f')

        if value1 == 'Poröser':
            value3 = column.number_input(f"Strömungswiderstand [Ns/m^4]", key=f"value3_{i}", format='%e')
            value4 = column.number_input(f"Porosität", key=f"value4_{i}", format='%0f', value=0.98)
            value5 = column.number_input(f"Tortuosität", key=f"value5_{i}", format='%0f', value=1.4)
            material_dict[key] = [value1, value2, value3, value4, value5]
        if value1 == 'Lochplatte':
            value6 = column.number_input(f"Lochdurchmesser [mm]", key=f"value6_{i}", format='%0f')
            value7 = column.number_input(f"Lochabstand [mm]", key=f"value7_{i}", format='%0f')
            material_dict[key] = [value1, value2, value6, value7]
        if value1 == 'Platte':
            value8 = column.number_input(f"Materialdichte [kg/m^3]", key=f"value8_{i}", format='%0f')
            value9 = column.number_input(f"Elastizitätsmodul [Pa]", key=f"value9_{i}", format='%e', value=4.1 * 10 ** 9)
            value10 = column.number_input(f"Nu", key=f"value10_{i}", format='%0f', value=0.3)
            value11 = column.number_input(f"Verlustfaktor", key=f"value11_{i}", format='%0f', value=0.1)
            material_dict[key] = [value1, value2, value8, value9, value10, value11]
        if value1 == 'Luft':
            material_dict[key] = [value1, value2]

    if st.button("OK"):
        st.write("Eingabe gespeichert!")
        st.write('  ' * 50000)
        st.write("Material parameter überprüfen:")
        for key, value in material_dict.items():
            st.write(f"{key}: {value}")

st.markdown('----')

################## Variable definition ##################
theta = theta * np.pi / 180
air_density = air_pressure / (287.058 * (air_temp + 273.15))
air_speed = 331.3 * np.sqrt(1 + (air_temp / 273.15))
viscosity = (1.458 * 10 ** (-6) * (air_temp + 273.15) ** (3 / 2)) / (air_temp + 273.15 + 110.4)
Z0 = air_speed * air_density
alphas = np.array([])

################## Computation ##################
try:
    for f in f_range_full:
        kx = 2 * np.pi * f / air_speed * np.sin(theta)
        T = []

        for l in range(num_materials):
            key = str("Material " + str(l + 1))
            L1 = material_dict[key][1] / 1000

            if material_dict[key][0] == 'Poröser':
                sigma = material_dict[key][2]
                phi = material_dict[key][3]
                alpha_inf = material_dict[key][4]
                T.append(models.Porous_Absorber_JAC(f, air_density, air_speed, L1, viscosity, sigma, air_pressure, phi,
                                                    alpha_inf, kx).get_T())
            if material_dict[key][0] == 'Lochplatte':
                d_hole = material_dict[key][2]
                a = material_dict[key][3]
                T.append(models.PerforatedPlate_Absorber(f, air_density, air_speed, L1, viscosity, d_hole, a).get_T())
            if material_dict[key][0] == 'Platte':
                density = material_dict[key][2]
                E = material_dict[key][3]
                nu = material_dict[key][4]
                eta = material_dict[key][5]
                T.append(
                    models.Plate_Absorber(f, air_density, air_speed, L1, viscosity, theta, density, E, nu, eta).get_T())
            if material_dict[key][0] == 'Luft':
                T.append(models.Air_Absorber(f, air_density, air_speed, L1, viscosity, kx).get_T())

        alphas = np.append(alphas, absorptioncoeff.AbsorptionCoeff(T, Z0, theta).abs_coeff())
except:
    pass

################## Output Section ##################
# Plotting
try:
    st.header('Plot :bar_chart:')
    if plot_type == 'Graph':
        fig1 = utils.plotly_go_line(x=f_range,
                                    y=alphas[f_min:f_max],
                                    x_label='Frequenz in [Hz]',
                                    y_label='Absorptionsgrad',
                                    title="Absorptionsgrad Plot")
        st.plotly_chart(fig1)

        # DF anzeigen
        col1, col2 = st.columns(2)
        col1.subheader('Daten :books:')
        df = pd.DataFrame({'Frequenz [Hz]': f_range, 'Absorptionsgrad [1]': alphas[0:len(f_range)]})
        st.dataframe(df, height=210)
        col2.subheader('Herunterladen :arrow_heading_down:')
        with col2:
            export = utils.create_df_export_button(
                df=df,
                title=f"Absorptionsgrad Berechnung",
                ts=None,
            )
    elif plot_type == 'Oktavbänder':
        fig1 = utils.plotly_freq_bands(x=f_range_full,
                                       y=alphas,
                                       x_label='Frequenz in [Hz]',
                                       y_label='Absorptionsgrad',
                                       title="Absorptionsgrad Oktavbänder",
                                       plot_type="oct")
        st.plotly_chart(fig1)

        # DF anzeigen
        col1, col2 = st.columns(2)
        col1.subheader('Daten :books:')
        df = pd.DataFrame({'Frequenz [Hz]': f_range_full, 'Absorptionsgrad [1]': alphas})
        st.dataframe(df, height=210)
        col2.subheader('Herunterladen :arrow_heading_down:')
        with col2:
            export = utils.create_df_export_button(
                df=df,
                title=f"Absorptionsgrad Berechnung",
                ts=None,
            )
    elif plot_type == 'Terzbänder':
        fig1 = utils.plotly_freq_bands(x=f_range_full,
                                       y=alphas,
                                       x_label='Frequenz in [Hz]',
                                       y_label='Absorptionsgrad',
                                       title="Absorptionsgrad Terzbänder",
                                       plot_type="third")
        st.plotly_chart(fig1)

        # DF anzeigen
        col1, col2 = st.columns(2)
        col1.subheader('Daten :books:')
        df = pd.DataFrame({'Frequenz [Hz]': f_range_full, 'Absorptionsgrad [1]': alphas})
        st.dataframe(df, height=210)
        col2.subheader('Herunterladen :arrow_heading_down:')
        with col2:
            export = utils.create_df_export_button(
                df=df,
                title=f"Absorptionsgrad Berechnung",
                ts=None,
            )
except:
    pass
