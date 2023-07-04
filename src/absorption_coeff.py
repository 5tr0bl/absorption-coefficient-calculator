import models
import numpy as np
import utils

class AbsorptionCoeff:

    # soll das hier auch in utils.py eurer Meinung nach?
    # absorber_layers sollte auch eine Liste von Modellen sein können

    def __init__(self, L1, L2, f_min, f_max, absorber_layers: models.AbsorberModelInterface, air_density, air_c, theta):
        self.alphas = np.array([])
        self.L1 = L1
        self.L2 = L2
        self.f_min = f_min
        self.f_max = f_max
        self.absorber_layers = absorber_layers
        self.air_density = air_density # luft_dichte
        self.air_c = air_c # Luft_c
        self.impedance = self.air_density * self.air_c
        self.theta = theta # Einfallswinkel

    def absorption_coeffs(self):
        # Initialize List of TMMs and alphas with the size of the freq range
        T_total = self.alphas = [None for _ in range(int(self.f_max) - int(self.f_min))]
        
        for layer in range(len(self.absorber_layers)):

            freq_index = 0 # list indexing beginning at 0, not at f_min

            for f in range(int(self.f_min), int(self.f_max)):
                try:
                    self.absorber_layers[layer].set_f(f) 
                    Z1 = self.absorber_layers[layer].get_Zp()
                    k1 = self.absorber_layers[layer].get_kp()
                    k0 = k2 = 2 * np.pi * f / self.air_c
                    
                    # Berechnung der TMM-Matrix
                    if not T_total[freq_index]:
                        T_total[freq_index] = utils.tmm(self.theta,k0,k1,self.L1,Z1,k2,self.L2,self.impedance)
                    else:
                        T_total[freq_index] = np.matmul(T_total[freq_index], utils.tmm(self.theta,k0,k1,self.L1,Z1,k2,self.L2,self.impedance))

                    # Absorptionsgrad berechnen -> gehört das in tmm-fkt?
                    R = (T_total[freq_index][0, 0] - T_total[freq_index][1, 0]*self.impedance) / (T_total[freq_index][0, 0] + T_total[freq_index][1, 0]*self.impedance) # Reflexionskoeffizient
                    alpha = 1 - np.abs(R) ** 2
                    self.alphas[freq_index] = alpha
                    freq_index += 1

                except (ZeroDivisionError, ValueError):
                    print('bitte richtige Zahl eingeben')

        return self.alphas