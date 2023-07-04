import numpy as np

class AbsorberModelInterface:
    # an interface upon which each absorber model will be implemented

    f = dichte = phi = alpha_unend = sigma = gamma = P0 =  viskosität_L = \
         thermisch_L = Pr = viskosität = 0.0
    
    def __init__(self, f, dichte, phi, alpha_unend, sigma, gamma, P0,  viskosität_L, \
         thermisch_L, Pr, viskosität):
        self.f = f
        self.dichte = dichte
        self.phi = phi
        self.alpha_unend = alpha_unend
        self.sigma = sigma
        self.gamma = gamma
        self.P0 = P0
        self.viskosität_L = viskosität_L
        self.thermisch_L = thermisch_L
        self.Pr = Pr
        self.viskosität = viskosität

    def calculate_aux_values(self):
        pass
    
    def get_density_p(self):
        pass
    
    def get_kp(self):
        pass

    def get_Zp(self):
        pass

    def set_f(self, new_f):
        self.f = new_f

# Helmholtz Resonator
# Acoustic Waves p.131
# Lossless model Z0=Zn=Zc, k0=kn=kc
# n = neck, c = cavity

class HelmholtzResonator(AbsorberModelInterface):
    def __init__(self, l_c, l_n):
        pass
    '''
    def __init__(self, f, dichte, phi, alpha_unend, sigma, gamma, P0,  viskosität_L,
                 thermisch_L, Pr, viskosität):
        super().__init__(f, dichte, phi, alpha_unend, sigma, gamma, P0, viskosität_L,
                         thermisch_L, Pr, viskosität)
    '''    
    def calculate_aux_values(self):
        pass
    
    def get_density_p(self):
        pass
    
    def get_kp(self):
        pass

    def get_Zp(self):
        pass

    def get_res_freq():
        pass

    def set_f(self, new_f):
        self.f = new_f

# Poröser Absorber
# Acoustic Waves p.128
# JAC-Modell
# Porosität des porösen Materials (Φ)
# Dichte der Flüssigkeit (normalerweise Luft) (ρ)
# Viskosität des Fluids (η)
# Strömungswiderstand des porösen Materials (σ)
# Frequenz der Schallwelle (f)
# Thermische Parameter eines porösen Materials (Pr, Le)
# Wärmekapazität des porösen Materials (Cp)
# Porengröße des porösen Materials (d)
class Porous_Absorber(AbsorberModelInterface):

    def __init__(self, f, dichte, phi, alpha_unend, sigma, gamma, P0,  viskosität_L,
                 thermisch_L, Pr, viskosität):
        super().__init__(f, dichte, phi, alpha_unend, sigma, gamma, P0, viskosität_L,
                         thermisch_L, Pr, viskosität)
        
        self.K0 = self.gamma * self.P0
        self.calculate_aux_values()
    
    def calculate_aux_values(self):
        try:
            self.omega = 2 * np.pi * self.f
            self.G1 = self.sigma * self.phi / (self.alpha_unend * self.dichte * self.omega)
            self.G2 = 4*((self.alpha_unend)**2) * self.dichte * self.viskosität * self.omega / ((self.sigma*self.phi*self.viskosität_L)**2)

            self.G1_dot = 8*self.viskosität / (self.dichte * self.Pr * ((self.thermisch_L)**2) * self.omega)

            self.G2_dot = self.dichte * self.Pr * ((self.thermisch_L)**2) * self.omega / (16*self.viskosität)

            self.dichte_p = self.dichte * self.alpha_unend * (1- 1j*self.G1*np.sqrt(1+1j*self.G2)) / self.phi

            self.K_p = self.K0*self.phi**(-1) / (self.gamma - (self.gamma-1)*((1- 1j*self.G1_dot*np.sqrt(1+ 1j*self.G2_dot))**-1))
        except ZeroDivisionError:
            print('bitte richtige Zahl eingeben')

    def get_density_p(self):
        return self.dichte_p
    
    def get_kp(self):
        return self.omega * np.sqrt(self.dichte_p/self.K_p)
    
    def get_Zp(self):
        return np.sqrt(self.dichte_p*self.K_p)
    
    def set_f(self, new_f):
        super().set_f(new_f)
        self.calculate_aux_values()


def JAC(f: float,
        dichte: float,
        phi: float, 
        alpha_unend: float,
        sigma: float,
        gamma: float,
        P0: float,
        viskosität_L: float,
        thermisch_L: float,
        Pr: float,
        viskosität: float
        ) -> float:

    K0 = gamma * P0
    omega = 2*np.pi*f
    G1 = sigma * phi / (alpha_unend * dichte * omega)

    G2 = 4*((alpha_unend)**2) * dichte * viskosität * omega / ((sigma*phi*viskosität_L)**2)

    G1_dot = 8*viskosität / (dichte * Pr * ((thermisch_L)**2) * omega)

    G2_dot = dichte * Pr * ((thermisch_L)**2) * omega / (16*viskosität)

    dichte_p = dichte * alpha_unend * (1- 1j*G1*np.sqrt(1+1j*G2)) / phi

    K_p = K0*phi**(-1) / (gamma - (gamma-1)*((1- 1j*G1_dot*np.sqrt(1+ 1j*G2_dot))**-1))

    kp = omega * np.sqrt(dichte_p/K_p)
    Zp = np.sqrt(dichte_p*K_p)

    return Zp, kp


