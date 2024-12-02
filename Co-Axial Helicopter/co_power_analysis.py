import math
import numpy as np
import matplotlib.pyplot as plt
from co_rotor_sizing import RotorSizing

class PowerAnalysis:

    def __init__(self):
        #conversion
        self.g = 9.80665
        self.pi = math.pi
        self.RPM_to_rad = (math.pi / 30)

        #input
        self.rotorsizing = RotorSizing()
        self.solidity = self.rotorsizing.maximum_solidity
        self.rho = 1.225
        self.omega = self.rotorsizing.omega
        self.rotor_radius = self.rotorsizing.rotor_radius
        self.mtow = self.rotorsizing.MTOW
        self.mtow_N = self.mtow * self.g
        self.k = 1.15
        self.k_dl = 1.04
        self.V_point = 13
        self.k_int = 1.28
        self.ROC_VCD = 0

        # basic gamma values
        self.gamma_CD = 0

        self.forward_flight()

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

        ## parasitic power calculation
        self.A_eq = 0.75
        self.P_par = 0.5*self.A_eq*self.rho*(self.V_point**3)
        self.D_par = self.P_par / self.V_point

        ## induced power calculation
        self.T = self.k_dl * self.mtow_N
        self.v_i_hov = math.sqrt((self.mtow_N)/(2*self.rho*self.pi*(self.rotor_radius**2)))
        self.V_bar = self.V_point / self.v_i_hov
        self.alpha = math.asin((self.D_par/self.mtow_N) + math.sin(math.radians(self.gamma_CD)))
        array = np.array([0.0,0.0,0.0,0.0,0.0])
        array[0] = 1
        array[1] = 2*self.V_bar*math.sin(self.alpha)
        array[2] = (self.V_bar**2)
        array[3] = 0
        array[4] = -1
        self.roots = np.roots(array)[3]
        self.v_i_bar = np.real(self.roots) #aanname
        self.v_i = self.v_i_bar*self.v_i_hov

        #compute induced power
        self.P_i = self.k_int * self.k * self.T * self.v_i

        ## total power
        self.P_total_level_loss = self.P_p + self.P_i + self.P_par

        # account for losses
        self.P_total_level = 1.045*self.P_total_level_loss

        ## compute climb/descent power
        #1.045 is a drag to fuselage kinda constant. 
        self.ROC_CD = self.V_point * math.sin(math.radians(self.gamma_CD))
        self.P_CD_loss = self.mtow_N * self.ROC_CD
        self.P_CD = self.P_CD_loss * 1.045

        ## total power (including climb/descent)
        self.P_total_CD = self.P_total_level + self.P_CD

        ## hover power
        self.P_hoge = self.k_int * self.k * self.T * self.v_i_hov + self.P_p_hov

        ## vertical climb/descent power
        self.P_VCD = self.P_hoge + self.ROC_VCD*self.mtow_N


    def iterate_design(self, new_mtow_N=None, new_V_point=None, new_solidity=None, new_gamma_CD=None, new_rho=None, new_ROC_VCD=None):
        if new_mtow_N:
            self.mtow_N = new_mtow_N
        if new_V_point:
            self.V_point = new_V_point
        if new_solidity:
            self.solidity = new_solidity
        if new_gamma_CD:
            self.gamma_CD = new_gamma_CD
        if new_rho:
            self.rho = new_rho
        if new_ROC_VCD:
            self.ROC_VCD = new_ROC_VCD


        self.forward_flight()

    def display_parameters(self):
        print("----------------")


if __name__ == '__main__':
    power = PowerAnalysis()
    power.display_parameters()

    #specify athmospheric condition
    power.iterate_design(new_rho=1.19011)

    ## Compute hover specifics
    # specify hover flight conditions
    power.iterate_design(new_ROC_VCD=0)
    print(f"HOGE power = {power.P_hoge} [W]")
    print(f"Vertical climb/descent power = {power.P_VCD} [W]")

    ## Compute forward flight specifics
    # specify forward flight condition (through gamma)
    power.iterate_design(new_gamma_CD=0)

    # initiate
    V = np.linspace(0.01, 85, 1000)
    AV = []
    P_p = []
    P_i = []
    P_par = []
    P_total_level = []
    P_CD = []
    P_total_CD = []

    for velocity in V:
        power.iterate_design(new_V_point=velocity)
        AV.append(velocity/(power.omega*power.rotor_radius))
        P_p.append(power.P_p)
        P_i.append(power.P_i)
        P_par.append(power.P_par)
        P_total_level.append(power.P_total_level)
        P_CD.append(power.P_CD)
        P_total_CD.append(power.P_total_CD)

    # compute velocity at minimum power required
    V_be = V[P_total_CD.index(min(P_total_CD))]

    # compute power required at Vbe
    power.iterate_design(new_V_point=V_be)

    print(f"Maximum endurance velocity = {V_be} [m/s] = {3.6 * V_be} [km/h]")
    print(f"Maximum endurance power = {power.P_total_CD} [W]")
    print(power.solidity,power.omega,power.rotor_radius)

    print("----------------")

    ## power plots
    print = True
    if print == True:
        plt.plot(AV, P_p, label="Profile drag power", linestyle='-', color='b')
        plt.plot(AV, P_i, label="Induced power", linestyle='--', color='c')
        plt.plot(AV, P_par, label="Parasitic power", linestyle='-.', color='r')
        plt.plot(AV, P_total_level, label="Total power (level flight)", linestyle='-', color='k')

        # add textual labels
        plt.text(AV[-1], P_p[-1], 'P_p', color='black', va='center', ha='left')
        plt.text(AV[-1], P_i[-1], 'P_i', color='black', va='center', ha='left')
        plt.text(AV[-1], P_par[-1], 'P_par', color='black', va='center', ha='left')
        plt.text(AV[-1], P_total_level[-1], 'P_total_level', color='black', va='center', ha='left')

        #climb descent
        plt.plot(AV, P_CD, label='Climb power ($\gamma$ = 9$^\circ$)', linestyle='-', color='m')
        plt.plot(AV, P_total_CD, linestyle=':', label='Total power (climbing flight)', color='k')

        #add textual labels
        plt.text(AV[-1], P_CD[-1], 'P_C', color='black', va='center', ha='left')
        plt.text(AV[-1], P_total_CD[-1], 'P_total_CD', color='black', va='center', ha='left')

        plt.axvline(x=0.4569090434, color='gray', linestyle='-', linewidth=1, label='Never Exceed AV')


        plt.title("Power Components vs. Advance Ratio")
        plt.legend(loc='upper left', fontsize='large')
        plt.grid(color='gray', linestyle=':', linewidth=0.5)
        plt.xlabel("Advance Ratio [-]")
        plt.ylabel("Power [W]")
        plt.show()




