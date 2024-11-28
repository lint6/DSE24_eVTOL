import math
import numpy as np
import matplotlib as plt

class PowerAnalysis:

    def __init__(self):
        #conversion
        self.g = 9.80665
        self.pi = math.pi
        self.RPM_to_rad = (math.pi / 30)


        #input
        self.solidity = 0.03
        self.rho = 1.225
        self.omega = 43.061
        self.rotor_radius = 4.7836157
        self.mtow = 718.89
        self.mtow_N = self.mtow * self.g
        self.k = 1.15
        self.k_dl = 1.04
        self.V_point = 13

        #advance ratio
        self.AV = self.V_point/(self.omega*self.rotor_radius)

        self.forward_flight()
        self.climbing_flight()


    def forward_flight(self):
        ## profile drag power calculation
        self.C_T = self.mtow_N/(self.rho * self.pi * (self.rotor_radius**2) * ((self.rotor_radius * self.omega)**2))
        self.C_L_bar = (6.6 * self.C_T)/self.solidity
        self.C_l_alpha = 5.73 #rad-1
        self.alpha_m = self.C_L_bar / self.C_l_alpha

        # estimate C_D_p_bar
        self.C_D_p_bar_1 = 0.0087 - 0.0216 * self.alpha_m + 0.4 * (self.alpha_m**2) #bailey
        self.C_D_p_bar_2 = 0.011 + 0.4 * (self.alpha_m**2) #marinescu
        self.C_D_p_bar_3 = 0.009 + 0.73 * (self.alpha_m**2) #talbot

        # pick the highest
        self.C_D_p_bar = max(self.C_D_p_bar_1,self.C_D_p_bar_2,self.C_D_p_bar_3)

        # compute hover profile power
        self.P_p_hov = ((self.solidity*self.C_D_p_bar)/8) * self.rho *((self.omega*self.rotor_radius)**3)*self.pi*(self.rotor_radius**2)

        # compute forward flight profile power
        self.P_p = self.P_p_hov*(1 + 4.65*(self.AV**2))

        ## induced power calculation
        self.T = self.k_dl * self.mtow_N
        self.v_i_hov = math.sqrt((self.mtow_N)/(2*self.rho*self.pi*(self.rotor_radius**2)))
        self.V_bar = self.V_point / self.v_i_hov
        self.v_i_bar = 1/self.V_bar
        self.v_i = self.v_i_hov * self.v_i_bar

        #compute induced power
        self.P_i = self.k * self.T * self.v_i

        ## parasite power calculation
        self.A_eq = 0.75
        self.P_par = 0.5*self.A_eq*self.rho*(self.V_point**3)

        ## total power calculation
        self.P_total_level_loss = self.P_p + self.P_i + self.P_par

        #account for losses
        self.P_total_level = 1.045*self.P_total_level_loss

    def climbing_flight(self):
        # compute climb power
        self.ROC = 2 #[m/s]
        self.P_C = self.mtow_N * self.ROC

        #total power
        self.P_total_climb = self.P_total_level + self.P_C



    def display_parameters(self):
        print(f"P_p = {self.P_p}")
        print(f"P_i = {self.P_i}")
        print(f"P_par = {self.P_par}")
        print(f"P_total_level = {self.P_total_level}")
        print(f"P_C = {self.P_C}")
        print(f"P_total_climb = {self.P_total_climb}")


if __name__ == '__main__':
    power = PowerAnalysis()
    power.display_parameters()
