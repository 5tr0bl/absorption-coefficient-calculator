import numpy as np
from scipy.special import jv


class AbsorberModelInterface:
    """Base class interface for all Absorber Models.

    Args:
        f (float): Frequency
        air_density (float): Density of air
        air speed (float): Speed of air
        L (float): Thickness of the layer
        viscosity (float): Viscosity of air
        omega (float): Angular frequency
    """

    f = air_density = omega = 0.0

    def __init__(self, f, air_density, air_speed, L1, viscosity):
        self.f = f
        self.air_density = air_density
        self.air_speed = air_speed
        self.L1 = L1
        self.omega = 2 * np.pi * self.f
        self.viscosity = viscosity

    def get_k(self):
        """Calculates the wave number. Different for each model, see source code for details.
        
        Returns:
            k (float): Wave number
        """
        pass

    def get_Z(self):
        """Calculates the surface impedance. Different for each model, see source code for details.

        Returns:
            Z (float): Surface impedance
        """
        pass

    def get_T(self):
        """Calculates the transfer matrix. Different for each model, see source code for details.

        Returns:
            T (float): Transfer Matrix of the absorber
        """
        pass


class Porous_Absorber_DB(AbsorberModelInterface):
    """Delany & Bazley Empirical Model for a porous absorber material.

    Args:
        f (float): Frequency
        air_density (float): Density of air
        air speed (float): Speed of air
        sigma (float): Flow resistivity of material
        L (float): Thickness of the layer
        viscosity (float): Viscosity of air
        kx (float): Wave number in x direction


    Returns:
        k (float): Wave number when calling get_k()
        Z (float): Surface impedance when calling get_Z()
        T (float): Transfer Matrix of the absorber when calling get_T()
    """

    def __init__(self, f, air_density, air_speed, L1, viscosity, sigma, kx):
        super().__init__(f, air_density, air_speed, L1, viscosity)

        self.sigma = sigma
        self.kx = kx
        self.X = (self.air_density * self.f) / self.sigma

    def get_k(self):
        k = self.omega / self.air_speed * (1 + 0.0978 * self.X ** (-0.7) - 1j * 0.189 * self.X ** (-0.595))
        return k

    def get_Z(self):
        Z = self.air_density * self.air_speed * (1 + 0.0571 * self.X ** (-0.754) - 1j * 0.087 * self.X ** (-0.732))
        return Z

    def get_T(self):
        k = self.get_k()
        Z = self.get_Z()
        k_z = np.sqrt(k ** 2 - self.kx ** 2)

        T = np.array([[np.cos(k_z * self.L1), 1j * Z * (k / k_z) * np.sin(k_z * self.L1)],
                      [(1j / Z) * (k_z / k) * np.sin(k_z * self.L1), np.cos(k_z * self.L1)]])
        return T


class Porous_Absorber_JAC(AbsorberModelInterface):
    """Johnson-Champoux-Allard Model for a porous absorber material.

    Args:
        f (float): Frequency
        air_density (float): Density of air
        air speed (float): Speed of air
        sigma (float): Flow resistivity of material
        L (float): Thickness of the layer
        viscosity (float): Viscosity of air

        air_pressure (float): Air pressure
        phi (float): Porosity
        alpha_inf (float): Tortuosity
        kx (float): Wave number in x direction

        gamma (float): Specific heat ratio
        kappa (float): Thermal conductivity
        cp (float): Specific heat capacity
        K0 (float): Bulk modulus
        delta_v (float): Viscous boundary layer thickness
        delta_h (float): Thermal boundary layer thickness
        Pr (float): Prandtl number


    Returns:
        k (float): Wave number
        Z (float): Surface impedance
        T (float): Transfer Matrix of the absorber
    """

    def __init__(self, f, air_density, air_speed, L1, viscosity, sigma, air_pressure, phi, alpha_inf, kx):
        super().__init__(f, air_density, air_speed, L1, viscosity)

        self.sigma = sigma
        self.air_pressure = air_pressure
        self.phi = phi
        self.alpha_inf = alpha_inf
        self.kx = kx
        self.gamma = 1.4
        self.K0 = self.gamma * self.air_pressure
        self.kappa = 0.0241
        self.cp = 1.01
        self.delta_v = np.sqrt(2 * self.viscosity / (self.air_density * self.omega))
        self.delta_h = np.sqrt(2 * self.kappa / (self.air_density * self.omega * self.cp))
        self.Pr = (self.delta_v / self.delta_h) ** 2
        self.viscosity_L = 1 / 1 * np.sqrt(8 * self.viscosity * self.alpha_inf / (self.phi * self.sigma))
        self.thermal_L = 2 * self.viscosity_L  # simpler formulations

        self.G1 = self.sigma * self.phi / (self.alpha_inf * self.air_density * self.omega)
        self.G2 = 4 * ((self.alpha_inf) ** 2) * self.air_density * self.viscosity * self.omega / ((self.sigma * self.phi * self.viscosity_L) ** 2)

        self.G1_dot = 8 * self.viscosity / (self.air_density * self.Pr * ((self.thermal_L) ** 2) * self.omega)
        self.G2_dot = self.air_density * self.Pr * ((self.thermal_L) ** 2) * self.omega / (16 * self.viscosity)
        self.density_p = self.air_density * self.alpha_inf * (1 - 1j * self.G1 * np.sqrt(1 + 1j * self.G2)) / self.phi
        self.Kp = self.K0 * self.phi ** (-1) / (self.gamma - (self.gamma - 1) * ((1 - 1j * self.G1_dot * np.sqrt(1 + 1j * self.G2_dot)) ** -1))

    def get_k(self):
        k = self.omega * np.sqrt(self.density_p / self.Kp)
        return k

    def get_Z(self):
        Z = np.sqrt(self.density_p * self.Kp)
        return Z

    def get_T(self):
        k = self.get_k()
        Z = self.get_Z()
        k_z = np.sqrt(k ** 2 - self.kx ** 2)

        T = np.array([[np.cos(k_z * self.L1), 1j * Z * (k / k_z) * np.sin(k_z * self.L1)],
                      [(1j / Z) * (k_z / k) * np.sin(k_z * self.L1), np.cos(k_z * self.L1)]])
        return T


class PerforatedPlate_Absorber(AbsorberModelInterface):
    """Maa´s Model for a rigid micro-perforated plate absorber material.

    Args:
        f (float): Frequency
        air_density (float): Density of air
        air speed (float): Speed of air
        L (float): Thickness of the layer
        viscosity (float): Viscosity of air

        d_hole (float): Diameter of hole
        a (float): Distance between holes

        phi (float): Porosity
        s (float): Ratio of the holes' diameter to the boundary layer thickness
        S (float): Area of the edge effects
        e (float): DONT KNOW THE NAME
        F_e (float): Fok function

    Returns:
        Z (float): Surface impedance
        k (float): Wave number
        T (float): Transfer Matrix of the absorber
    """

    def __init__(self, f, air_density, air_speed, L1, viscosity, d_hole, a):
        super().__init__(f, air_density, air_speed, L1, viscosity)

        self.d_hole = d_hole/1000
        self.a = a/1000

    def get_k(self):
        k = self.omega / self.air_speed
        return k

    def get_Z(self):
        self.phi = (np.pi / 4) * (self.d_hole / self.a) ** 2
        self.e = 1.1284 * np.sqrt(self.phi)
        self.s = self.d_hole * np.sqrt(self.air_density * self.omega / 4 / self.viscosity)
        self.F_e = (1 - 1.4092 * self.e + 0.33818 * (self.e ** 3) + 0.06793 *
                    (self.e ** 5) - 0.02287 * (self.e ** 6) + 0.03015 *
                    (self.e ** 7) - 0.01641 * (self.e ** 8)) ** (-1)
        J_0 = jv(0, self.s * np.sqrt(-1j))
        J_1 = jv(1, self.s * np.sqrt(-1j))
        Z = ((np.sqrt(2 * self.air_density * self.omega * self.viscosity) / 2 * self.phi) +
             (1j * (self.omega * self.air_density / self.phi)) * (0.85 * self.d_hole / self.F_e +
                                                self.L1 * (1 - 2 * J_1 / (self.s * np.sqrt(-1j)) / J_0) ** (-1)))
        return Z

    def get_T(self):
        T = np.array([[1, self.get_Z()],
                      [0, 1]])
        return T


class Air_Absorber(AbsorberModelInterface):
    """Air Model

    Args:
        f (float): Frequency
        air_density (float): Density of air
        air speed (float): Speed of air
        L (float): Thickness of the layer
        viscosity (float): Viscosity of air

        kx (float): Wave number in x direction

    Returns:
        Z (float): Surface impedance
        k (float): Wave number
        T (float): Transfer Matrix of the absorber
    """

    def __init__(self, f, air_density, air_speed, L1, viscosity, kx):
        super().__init__(f, air_density, air_speed, L1, viscosity)

        self.kx = kx

    def get_k(self):
        k = 2 * np.pi * self.f / self.air_speed
        return k

    def get_Z(self):
        Z = self.air_density * self.air_speed
        return Z

    def get_T(self):
        k = self.get_k()
        Z = self.get_Z()
        k_z = np.sqrt(k ** 2 - self.kx ** 2)
        T = np.array([[np.cos(k_z * self.L1), 1j * Z * (k / k_z) * np.sin(k_z * self.L1)],
                      [(1j / Z) * (k_z / k) * np.sin(k_z * self.L1), np.cos(k_z * self.L1)]])
        return T

class Plate_Absorber(AbsorberModelInterface):
    """Infinite Elastic Vibrating Wall Model for a plate absorber material.

        Args:
            f (float): Frequency
            air_density (float): Density of air
            air speed (float): Speed of air
            L (float): Thickness of the layer
            viscosity (float): Viscosity of air

            theta (float): Angle of incidence
            density (float): Density of the plate
            E (float): Young's modulus
            nu (float): Poisson's ratio
            eta (float): Loss factor


        Returns:
            fc (float): Critical frequency
            Z (float): Surface impedance
            T (float): Transfer Matrix of the absorber
        """

    def __init__(self, f, air_density, air_speed, L1, viscosity, theta, density, E, nu, eta):
        super().__init__(f, air_density, air_speed, L1, viscosity)

        self.omega = 2 * np.pi * self.f
        self.theta = theta
        self.density = density
        self.E = E
        self.nu = nu
        self.eta = eta

        self.m_dot = self.density * self.L1
        self.D = self.E * self.L1 ** 3 / (12 * (1 - self.nu ** 2))

    def get_fc(self):
        fc = (self.air_speed ** 2) / (2 * np.pi) * np.sqrt(self.m_dot / self.D)
        return fc

    def get_Z(self):
        fc = self.get_fc()
        Z = 1j * self.m_dot * self.omega * (1 - ((self.f/fc)**2) * (1 + 1j*self.eta) * np.sin(self.theta)**4)
        return Z

    def get_T(self):
        Z = self.get_Z()
        T = np.array([[1, Z],
                      [0, 1]])
        return T
