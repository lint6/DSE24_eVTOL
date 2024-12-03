import math
import numpy as np
import matplotlib.pyplot as plt

class RotorSizing:

    def __init__(self, MTOW=718.89, n_blades=4, n_rotor = 4, bank_angle = 30, cto_fl = 0.12, cto_turn = 0.15, cto_turb = 0.17, d_fact = 0.05):

        # conversions  
        self.celsius_to_kelvin = 273.15     # add
        self.RPM_to_rad = math.pi / 30      # multiply
        self.deg_to_rad = math.pi / 180     # multiply
        self.ft_to_m = 0.3048               # multiply

        # ambient constants 
        self.g = 9.80665        # m / s^2
        self.rho = 1.225        # kg/m^3
        self.temperature = 15   # celsius 
        self.speed_of_sound = math.sqrt(1.4 * 287 * (self.temperature + self.celsius_to_kelvin))

        # physical constants
        self.MTOW = MTOW 
        self.n_blades = n_blades 
        self.max_tip_mach = 0.85
        self.N_rotors = n_rotor
        self.c_t_o_fl = cto_fl #user defined depending on advance ratio
        self.c_t_o_turn = cto_turn # user defined
        self.c_t_o_turb = cto_turb # user defined
        self.roll_angle = bank_angle                                        # deg
        self.n_z_turn = 1 / np.cos(self.roll_angle * self.deg_to_rad)    # load factor in turn
        self.lift_slope = 5.73      # 1 / rad (from NACA0012)
        self.gust_velocity = 30 * self.ft_to_m # from FAA
        self.coaxial = 1 # 1 if not coaxial, 2 if coaxial

        download_factor = d_fact
        self.fuselage_download = download_factor * self.MTOW * self.g                 # N
        self.k_dl = 1 + (self.fuselage_download / (self.MTOW * self.g))    # --

        # initialize the parameters 
        self.update_parameters()

    def update_parameters(self):
        # derived parameters 
        # Disc loading calcualtion, we assumed the disc loading relation, To be further discussed 
        #check the source of this 
        #for co-axial
        #self.disc_loading = 2.28 * ((self.MTOW)**(1/3) - 2.34)                      # kg/m^2

        # for quad (assumed):
        self.disc_loading = 34.2   # kg/m^2 ( 7 lbs/ft2)
        self.rotor_radius = np.sqrt((self.MTOW / self.N_rotors) / (np.pi * self.disc_loading))    # m 
        self.rotor_diameter = 2 * self.rotor_radius                                 # m
        # marilena statistics for helicopters. quad assumed 550 ft/s (NASA paper), to be discussed
        #self.tip_speed = 140 * (self.rotor_diameter)**0.171                         # m / s
        self.tip_speed = 550 * self.ft_to_m

        #self.RPM = 2673 / ((self.rotor_diameter)**0.892)                            # rpm
        #self.omega = self.RPM_to_rad * self.RPM                                     # rad / s

        self.omega = self.tip_speed / self.rotor_radius # rad/s
        self.RPM = self.omega * (60/2*np.pi) # rpm

        # flight performance 
        self.max_forward_velocity = (self.max_tip_mach * self.speed_of_sound) - self.tip_speed  # m / s; choose max mach, calculate max velocity
        self.never_exceed_velocity = 1.1 * self.max_forward_velocity                            # m / s
        self.mu_ne = self.never_exceed_velocity / self.tip_speed                                # --; advance ratio

        # forward flight solidity 
        self.T_forward_flight = self.k_dl * self.MTOW * self.g             # N
        self.C_T_forward_flight = self.T_forward_flight / (self.N_rotors * (self.rho * math.pi * self.rotor_radius**2 * (self.rotor_radius * self.omega)**2)) # --
        self.solidity_forward_flight = self.C_T_forward_flight / self.c_t_o_fl                    # --

        # turn solidity 
        self.T_turn = self.n_z_turn * self.k_dl * self.MTOW * self.g     # N 
        self.C_T_turn = self.T_turn / (self.N_rotors*(self.rho * math.pi * (self.rotor_radius**2) * (self.rotor_radius * self.omega)**2))   # -- 
        self.solidity_turn = self.C_T_turn / self.c_t_o_turn

        # turbulent solidity
        self.delta_n_turbulent = (0.25 * self.lift_slope * (self.gust_velocity/(self.rotor_radius * self.omega))) / self.c_t_o_turb
        self.n_z_turbulent = self.delta_n_turbulent + 2 # 2g pull up
        self.T_turbulence = self.n_z_turbulent * self.k_dl * self.MTOW * self.g
        self.C_T_turbulence = self.T_turbulence / (self.N_rotors * (self.rho * math.pi * (self.rotor_radius**2) * (self.rotor_radius * self.omega)**2))
        self.solidity_turbulent = self.C_T_turbulence / self.c_t_o_turb

        # chord and aspect ratio 
        self.maximum_solidity = np.max([self.solidity_forward_flight, self.solidity_turn, self.solidity_turbulent])
        self.chord = (self.maximum_solidity * np.pi * self.rotor_radius) / (self.coaxial * self.n_blades)
        self.aspect_ratio = self.rotor_radius**2 / (self.rotor_radius * self.chord)


    def display_parameters(self):
        print(f"MTOW: {self.MTOW} kg")
        print(f"Rotor Radius: {self.rotor_radius:.2f} m")
        print(f"Rotor Diameter: {self.rotor_diameter:.2f} m")
        print(f"Tip Speed: {self.tip_speed:.2f} m/s")
        print(f"RPM: {self.RPM:.2f}")
        print(f"Maximum Solidity: {self.maximum_solidity:.3f}")
        print(f"Blade Chord: {self.chord:.3f} m")
        print(f"Aspect Ratio: {self.aspect_ratio:.2f}")
        print(f"Max Forward Velocity: {self.max_forward_velocity:.2f} m/s")
        print(f"Never Exceed Velocity: {self.never_exceed_velocity:.2f} m/s")

    def iterate_design(self, new_MTOW=None, new_n_blades=None):
        if new_MTOW:
            self.MTOW = new_MTOW 
        if new_n_blades:
            self.n_blades = new_n_blades 
        self.update_parameters()

    def visual_blade_vs_aspect_ratio(self):
        blade_numbers = [2, 3, 4, 5, 6, 7, 8]
        aspect_ratios = []
        for blades in blade_numbers:
            self.iterate_design(new_n_blades=blades)
            aspect_ratios.append(self.aspect_ratio)
        plt.figure(figsize=(10, 6))
        plt.plot(blade_numbers, aspect_ratios, marker='o', label='Aspect Ratio')
        plt.title('Effect of Blade Number on Aspect Ratio')
        plt.xlabel('Number of Blades')
        plt.ylabel('Aspect Ratio')

        plt.axhspan(14, 20, color='green', alpha=0.2, label='Acceptable Range')
        plt.axhspan(0, 14, color = 'red', alpha=0.2, label = 'Unacceptable Range')
        plt.axhspan(20, 100, color = 'red', alpha = 0.2)

        plt.grid(True)
        plt.legend()
        plt.show()



class PowerAnalysis:

    def __init__(self, FM = 0.7, k_int = 1, A_eq = 0.418):
        #conversion
        self.g = 9.80665
        self.pi = np.pi
        self.RPM_to_rad = (np.pi / 30)

        #input
        self.rotorsizing = RotorSizing()
        self.FM = FM
        self.solidity = self.rotorsizing.maximum_solidity
        self.rho = self.rotorsizing.rho
        self.omega = self.rotorsizing.omega
        self.rotor_radius = self.rotorsizing.rotor_radius
        self.MTOW = self.rotorsizing.MTOW
        self.MTOW_N = self.MTOW * self.g
        self.k = 1.15       # inflow correction, assumption: should be between 1.1 and 1.2
        self.k_dl = self.rotorsizing.k_dl # fuselage downlaod factor
        self.V_point = 13 # cruise speed, assumption for now
        self.k_int = k_int # 1.28 for co axial
        self.ROC_VCD = 10
        self.C_l_alpha = self.rotorsizing.lift_slope # 1/rad NACA0012
        self.A_eq = A_eq # equivalent flat plate area

        # basic gamma values
        self.gamma_CD = 0

        self.forward_flight()

    def forward_flight(self):
        self.AV = self.V_point/(self.omega*self.rotor_radius) # advance ratio

        ## profile drag power calculation
        # thurst coefficient
        self.C_T = self.rotorsizing.T_forward_flight / (self.rho * (self.pi * (self.rotor_radius**2)) * self.rotorsizing.N_rotors * ((self.rotor_radius * self.omega)**2))
        # avg coefficient of lift
        self.C_L_bar = (6.6 * self.C_T)/self.solidity
        # AoA for max lift
        self.alpha_m = self.C_L_bar / self.C_l_alpha

        # estimate C_D_p_bar
        self.C_D_p_bar_1 = 0.0087 - 0.0216 * self.alpha_m + 0.4 * (self.alpha_m**2) #bailey
        self.C_D_p_bar_2 = 0.011 + 0.4 * (self.alpha_m**2) #marinescu
        self.C_D_p_bar_3 = 0.009 + 0.73 * (self.alpha_m**2) #talbot

        # pick the highest
        self.C_D_p_bar = max(self.C_D_p_bar_1,self.C_D_p_bar_2,self.C_D_p_bar_3)

        # compute hover profile power
        self.P_p_hov = ((self.solidity*self.C_D_p_bar)/8) * self.rho *((self.omega*self.rotor_radius)**3)*(self.pi*(self.rotor_radius**2)) * self.rotorsizing.N_rotors

        # compute forward flight profile power
        self.P_p = self.P_p_hov*(1 + 4.65*(self.AV**2))

        ## parasitic power calculation
        self.P_par = 0.5*self.A_eq*self.rho*(self.V_point**3)
        #parasitic drag
        self.D_par = self.P_par / self.V_point

        ## induced power calculation
        self.T = self.rotorsizing.T_forward_flight
        self.v_i_hov = math.sqrt((self.MTOW_N / self.rotorsizing.N_rotors)/(2*self.rho*self.pi*(self.rotor_radius**2)))
        self.V_bar = self.V_point / self.v_i_hov
        self.alpha = math.asin((self.D_par/self.MTOW_N) + math.sin(math.radians(self.gamma_CD)))
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
        self.P_CD_loss = self.MTOW_N * self.ROC_CD
        self.P_CD = self.P_CD_loss * 1.045

        ## total power (including climb/descent)
        self.P_total_CD = self.P_total_level + self.P_CD

        ## hover power; what is this?? changed it to momentum theory (james wang adn stepiniewski)
        #self.P_hoge = self.k_int * self.k * self.T * self.v_i_hov + self.P_p_hov

        self.P_hoge = ((self.MTOW * self.g) ** (3/2)) / (np.sqrt(2 * self.rho *(self.pi*(self.rotor_radius**2)) * self.rotorsizing.N_rotors) * self.FM)

        ## vertical climb/descent power
        self.P_VCD = self.P_hoge + self.ROC_VCD*self.MTOW_N


    def iterate_design(self, new_MTOW_N=None, new_V_point=None, new_solidity=None, new_gamma_CD=None, new_rho=None, new_ROC_VCD=None):
        if new_MTOW_N:
            self.MTOW_N = new_MTOW_N
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
    rotor = RotorSizing()
    print('----------------------------------------')
    print(f'Number of Blades = {rotor.n_blades} and MTOW = {rotor.MTOW}kg')
    print('----------------------------------------')
    rotor.display_parameters()
    ### visualizations ###
    rotor.visual_blade_vs_aspect_ratio()



if __name__ == '__main__':
    power = PowerAnalysis()
    power.display_parameters()

    #specify athmospheric condition
    power.iterate_design(new_rho=1.19011) # rho at cruise altitude

    ## Compute hover specifics
    # specify hover flight conditions
    power.iterate_design(new_ROC_VCD=0)
    print(f"HOGE power = {power.P_hoge / 1000:.2f} [kW]")
    print(f"Vertical climb/descent power = {power.P_VCD / 1000:.2f} [kW]")

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
    P_hoge = []

    for velocity in V:
        power.iterate_design(new_V_point=velocity)
        AV.append(velocity/(power.omega*power.rotor_radius))
        P_p.append(power.P_p / 1000)
        P_i.append(power.P_i / 1000)
        P_par.append(power.P_par / 1000)
        P_total_level.append(power.P_total_level / 1000)
        P_CD.append(power.P_CD / 1000)
        P_total_CD.append(power.P_total_CD / 1000)
        P_hoge.append(power.P_hoge / 1000)

    # compute velocity at minimum power required
    V_be = V[P_total_CD.index(min(P_total_CD))]

    # compute power required at Vbe
    power.iterate_design(new_V_point=V_be)

    print(f"Maximum endurance velocity = {V_be:.2f} [m/s] = {3.6 * V_be:.2f} [km/h]")
    print(f"Maximum endurance power = {power.P_total_level / 1000 :.2f} [kW]")
    print(f'Maximum solidity = {power.solidity:.3f}')
    print(f'Rotational Velocity = {power.omega:.2f} [rad/s]')
    print(f'Rotor Radius = {power.rotor_radius:.2f} [m]')

    print("----------------")

    ## power plots
    print = True
    if print == True:
        plt.plot(V, P_p, label="Profile drag power", linestyle='-', color='b')
        plt.plot(V, P_i, label="Induced power", linestyle='--', color='c')
        plt.plot(V, P_par, label="Parasitic power", linestyle='-.', color='r')
        plt.plot(V, P_total_level, label="Total power (level flight)", linestyle='-', color='k')

        # plot hover power
        plt.plot(V, P_hoge, label="Power HOGE", linestyle=':', color='purple')

        # add textual labels
        plt.text(V[-1], P_p[-1], 'P_p', color='black', va='center', ha='left')
        plt.text(V[-1], P_i[-1], 'P_i', color='black', va='center', ha='left')
        plt.text(V[-1], P_par[-1], 'P_par', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_level[-1], 'P_total_level', color='black', va='center', ha='left')
        plt.text(V[-1], P_hoge[-1], 'P_hoge', color='black', va='center', ha='left')

        #climb descent
        plt.plot(V, P_CD, label='Climb power ($\gamma$ = 9$^\circ$)', linestyle='-', color='m')
        plt.plot(V, P_total_CD, linestyle=':', label='Total power (climbing flight)', color='k')

        #add textual labels
        plt.text(V[-1], P_CD[-1], 'P_C', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_CD[-1], 'P_total_CD', color='black', va='center', ha='left')

        plt.axvline(x=0.4569090434, color='gray', linestyle='-', linewidth=1, label='Never Exceed AV')


        plt.title("Power Components vs. Advance Ratio")
        plt.legend(loc='upper left', fontsize='large')
        plt.grid(color='gray', linestyle=':', linewidth=0.5)
        plt.xlabel("Velocity [m/s]")
        plt.ylabel("Power [kW]")
        plt.show()