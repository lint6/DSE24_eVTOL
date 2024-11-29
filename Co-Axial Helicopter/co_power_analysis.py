import math
import numpy as np
import matplotlib.pyplot as plt

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

        # basic gamma values
        self.gamma_C = 0
        self.gamma_D = 0

        #advance ratio
        self.forward_flight()
        self.climbing_flight()
        self.descending_flight()
        self.iterate_design()


    def forward_flight(self):
        self.AV = self.V_point/(self.omega*self.rotor_radius)

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
        self.ROC_C = self.V_point * math.sin(math.radians(self.gamma_C))
        self.P_C = self.mtow_N * self.ROC_C

        #total power
        self.P_total_climb = self.P_total_level + self.P_C

    def descending_flight(self):
        # compute climb power
        self.ROC_D = self.V_point * math.sin(math.radians(self.gamma_D))
        self.P_D = self.mtow_N * self.ROC_D

        #total power
        self.P_total_descent = self.P_total_level + self.P_D

    def iterate_design(self, new_mtow_N=None, new_V_point=None, new_solidity=None, new_gamma_C=None, new_gamma_D=None):
        if new_mtow_N:
            self.mtow_N = new_mtow_N
        if new_V_point:
            self.V_point = new_V_point
        if new_solidity:
            self.solidity = new_solidity
        if new_gamma_C:
            self.gamma_C = new_gamma_C
        if new_gamma_D:
            self.gamma_D = new_gamma_D

        self.forward_flight()
        self.climbing_flight()
        self.descending_flight()


    def display_parameters(self):
        # print(f"P_p = {self.P_p}")
        # print(f"P_i = {self.P_i}")
        # print(f"P_par = {self.P_par}")
        print(f"P_total_level = {self.P_total_level}")
        print(f"P_C = {self.P_C}")
        # print(f"P_total_climb = {self.P_total_climb}")
        print(f"P_D = {self.P_D}")


if __name__ == '__main__':
    power = PowerAnalysis()
    #power.display_parameters()

    power.iterate_design(new_gamma_C=0)
    power.iterate_design(new_gamma_D=-9)

    V = np.linspace(5,100, 300)
    P_p = []
    P_i = []
    P_par = []
    P_total_level = []
    P_C = []
    P_total_climb = []
    P_D = []
    P_total_descent = []

    for velocity in V:
        power.iterate_design(new_V_point=velocity)
        P_p.append(power.P_p)
        P_i.append(power.P_i)
        P_par.append(power.P_par)
        P_total_level.append(power.P_total_level)
        P_C.append(power.P_C)
        P_total_climb.append(power.P_total_climb)
        P_D.append(power.P_D)
        P_total_descent.append(power.P_total_descent)

    climb = False
    descent = True


    plt.plot(V, P_p, linestyle='-', color='b')
    plt.plot(V, P_i, linestyle='--', color='c')
    plt.plot(V, P_par, linestyle='-.', color='r')
    plt.plot(V, P_total_level, linestyle='-', color='k')

    # add textual labels
    plt.text(V[-1], P_p[-1], 'P_p', color='black', va='center', ha='left')
    plt.text(V[-1], P_i[-1], 'P_i', color='black', va='center', ha='left')
    plt.text(V[-1], P_par[-1], 'P_par', color='black', va='center', ha='left')
    plt.text(V[-1], P_total_level[-1], 'P_total_level', color='black', va='center', ha='left')

    if climb:
        plt.plot(V, P_C, linestyle='-', color='m')
        plt.plot(V, P_total_climb, linestyle=':', color='k')

        plt.text(V[-1], P_C[-1], 'P_C', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_climb[-1], 'P_total_climb', color='black', va='center', ha='left')
    if descent:
        plt.plot(V, P_D, linestyle='-', color='g')
        plt.plot(V, P_total_descent, linestyle='--', color='k')

        plt.text(V[-1], P_D[-1], 'P_D', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_descent[-1], 'P_total_descent', color='black', va='center', ha='left')

    plt.title("Power Components vs. Velocity")
    plt.xlabel("Velocity (V) [m/s]")
    plt.ylabel("Power [W]")

    plt.show()

