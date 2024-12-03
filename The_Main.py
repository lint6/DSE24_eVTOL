import math
import numpy as np
import matplotlib.pyplot as plt

class RotorSizing:

    def __init__(self, MTOW=718.89, n_blades=4, n_rotor = 4, DL = 34.2, bank_angle = 30, cto_fl = 0.12, cto_turn = 0.15, cto_turb = 0.17, co_ax = 1, d_fact = 0.05, max_v = 50, k_int = 1, A_eq = 0.48, FM = 0.7):

        # conversions  
        self.celsius_to_kelvin = 273.15     # addition
        self.RPM_to_rad = np.pi / 30      # multiply
        self.deg_to_rad = np.pi / 180     # multiply
        self.ft_to_m = 0.3048               # multiply

        # ambient constants 
        self.g = 9.80665        # m / s^2
        self.rho = 1.225        # kg/m^3
        self.temperature = 15   # celsius 
        self.speed_of_sound = np.sqrt(1.4 * 287 * (self.temperature + self.celsius_to_kelvin))

        # physical constants
        self.MTOW = MTOW 
        self.n_blades = n_blades 
        self.max_tip_mach = 0.85
        self.N_rotors = n_rotor
        self.disc_loading = DL   # kg/m^2 
        self.c_t_o_fl = cto_fl #user defined depending on advance ratio
        self.c_t_o_turn = cto_turn # user defined
        self.c_t_o_turb = cto_turb # user defined
        self.roll_angle = bank_angle                                        # deg
        self.n_z_turn = 1 / np.cos(self.roll_angle * self.deg_to_rad)    # load factor in turn
        self.lift_slope = 5.73      # 1 / rad (from NACA0012)
        self.gust_velocity = 30 * self.ft_to_m # from FAA
        self.coaxial = co_ax # 1 if not coaxial, 2 if coaxial
        self.V_max = max_v
        self.V_ne = 1.1 * self.V_max

        download_factor = d_fact
        self.fuselage_download = download_factor * self.MTOW * self.g                 # N
        self.k_dl = 1 + (self.fuselage_download / (self.MTOW * self.g))    # --

        self.k_int = k_int
        self.A_eq = A_eq
        self.FM = FM

        # initialize the parameters 
        self.update_parameters()

    def update_parameters(self):
        # derived parameters 
        # Disc loading calcualtion, we assumed the disc loading relation, To be further discussed 
        #check the source of this 
        #for co-axial
        #self.disc_loading = 2.28 * ((self.MTOW)**(1/3) - 2.34)                      # kg/m^2

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

        #print(f'User chosen max velocity = {self.V_max * 3.6 :.2f} [km/h]')
        #print(f'Corresponding never exceed advance ratio = {self.V_ne / self.tip_speed} [-]')


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

        #plt.axhspan(14, 20, color='green', alpha=0.2, label='Acceptable Range')
        #plt.axhspan(0, 14, color = 'red', alpha=0.2, label = 'Unacceptable Range')
        #plt.axhspan(20, 100, color = 'red', alpha = 0.2)

        plt.grid(True)
        plt.legend()
        plt.show()



class PowerAnalysis:

    def __init__(self, rotorsizing):
        #conversion
        self.g = 9.80665
        self.pi = np.pi
        self.RPM_to_rad = (np.pi / 30)

        #input
        self.rotorsizing = rotorsizing if rotorsizing else RotorSizing()  # Use the provided object or create a default one
        self.FM = self.rotorsizing.FM
        self.solidity = self.rotorsizing.maximum_solidity
        self.rho = self.rotorsizing.rho
        self.omega = self.rotorsizing.omega
        self.rotor_radius = self.rotorsizing.rotor_radius
        self.MTOW = self.rotorsizing.MTOW
        self.MTOW_N = self.MTOW * self.g
        self.k = 1.15       # inflow correction, assumption: should be between 1.1 and 1.2
        self.k_dl = self.rotorsizing.k_dl # fuselage downlaod factor
        self.V_point = 13 # cruise speed, assumption for now
        self.k_int = self.rotorsizing.k_int # 1.28 for co axial
        self.ROC_VCD = 10
        self.C_l_alpha = self.rotorsizing.lift_slope # 1/rad NACA0012
        self.A_eq = self.rotorsizing.A_eq # equivalent flat plate area

        # basic gamma values
        self.gamma_CD = 9 # deg
        self.ROC_CD_steep = -7.6 # m/s

        self.hige_factor = 1  #power in hige / power hoge
        self.min_power = None
        self.min_power_velocity = None
        self.min_power_CD = None
        self.min_power_velocity_CD = None
        self.gamma_CD_steep = None

        self.P = {
            'HIGE1': None,
            'V_climb': None,
            'HOGE1': None,
            'Climb1': None,
            'Cruise1': None,
            'Descent1': None,
            'HOGE2': None,
            'Loiter': None,
            'Climb2': None,
            'Cruise2': None,
            'Descent2': None,
            'HOGE3': None,
            'V_Descent': None,
            'HIGE2': None
        }

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
        self.v_i_hov = np.sqrt((self.MTOW_N / self.rotorsizing.N_rotors)/(2*self.rho*self.pi*(self.rotor_radius**2)))
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

        #steep descent
        self.P_CD_loss_steep = self.MTOW_N * self.ROC_CD_steep
        self.P_CD_steep = self.P_CD_loss_steep * 1.045
        ## total power steep descent (including climb/descent)
        self.P_total_CD_steep = self.P_total_level + self.P_CD_steep

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


    def plot_power_components(self):
        V = np.linspace(0.01, 85, 1000)
        AV = []
        P_p = []
        P_i = []
        P_par = []
        P_total_level = []
        P_CD = []
        P_total_CD = []
        P_hoge = []
        P_total_CD_steep = []

        for velocity in V:
            self.iterate_design(new_V_point=velocity)
            AV.append(velocity / (self.omega * self.rotor_radius))
            P_p.append(self.P_p / 1000)
            P_i.append(self.P_i / 1000)
            P_par.append(self.P_par / 1000)
            P_total_level.append(self.P_total_level / 1000)
            P_CD.append(self.P_CD / 1000)
            P_total_CD.append(self.P_total_CD / 1000)
            P_hoge.append(self.P_hoge / 1000)
            P_total_CD_steep.append(self.P_total_CD_steep / 1000)

        # Identify velocity corresponding to minimum flight-level power
        self.min_power = min(P_total_level)
        self.min_power_velocity = V[P_total_level.index(self.min_power)]
        print(f"The velocity corresponding to the minimum flight-level power ({self.min_power:.2f} kW) is {self.min_power_velocity:.2f} m/s")

        #steep deescent angle
        #self.gamma_CD_steep = math.degrees(math.asin(self.ROC_CD_steep / self.min_power_velocity_CD))

        self.min_power_CD = min(P_total_CD)
        self.min_power_velocity_CD = V[P_total_CD.index(self.min_power_CD)]


        print(f"The velocity corresponding to the minimum climb power ({self.min_power_CD:.2f} kW) is {self.min_power_velocity_CD:.2f} m/s")
        plt.scatter(self.min_power_velocity_CD, self.min_power_CD, color='green', label='CD Min. Power', zorder=5)
        plt.annotate(f'({self.min_power_velocity_CD:.2f} m/s, {self.min_power_CD:.2f} kW)',
                 (self.min_power_velocity_CD, self.min_power_CD),
                 textcoords="offset points",
                 xytext=(-30, 10),
                 ha='center',
                 fontsize=9,
                 color='green')

        #plt.plot(V, P_p, label="Profile drag power", linestyle='-', color='b')
        #plt.plot(V, P_i, label="Induced power", linestyle='--', color='c')
        #plt.plot(V, P_par, label="Parasitic power", linestyle='-.', color='r')
        plt.plot(V, P_total_level, label="Total power (level flight)", linestyle='-', color='k')

        plt.plot(V, P_hoge, label="Power HOGE", linestyle=':', color='purple')

        # plotting min power point
        plt.scatter(self.min_power_velocity, self.min_power, color='red', label='Level Flight Min. Power', zorder=5)
        plt.annotate(f'({self.min_power_velocity:.2f} m/s, {self.min_power:.2f} kW)',
                 (self.min_power_velocity, self.min_power),
                 textcoords="offset points",
                 xytext=(-30, 10),
                 ha='center',
                 fontsize=9,
                 color='red')


        #plt.text(V[-1], P_p[-1], 'P_p', color='black', va='center', ha='left')
        #plt.text(V[-1], P_i[-1], 'P_i', color='black', va='center', ha='left')
        #plt.text(V[-1], P_par[-1], 'P_par', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_level[-1], 'P_total_level', color='black', va='center', ha='left')
        plt.text(V[-1], P_hoge[-1], 'P_hoge', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_CD_steep[-1], 'P_steep_descent', color='black', va='center', ha='left')

        #plt.plot(V, P_CD, label=f'Climb power ($\gamma$ = {self.gamma_CD}$^\circ$)', linestyle='-', color='m')
        plt.plot(V, P_total_CD, linestyle=':', label='Total power (climbing flight)', color='k')

        plt.plot(V, P_total_CD_steep, linestyle=':', label='Steep descent power', color='green')

        #plt.text(V[-1], P_CD[-1], 'P_C', color='black', va='center', ha='left')
        plt.text(V[-1], P_total_CD[-1], 'P_total_CD', color='black', va='center', ha='left')

        #plt.axvline(x=0.4569090434, color='gray', linestyle='-', linewidth=1, label='Never Exceed AV')

        plt.title("Power Components vs. Advance Ratio")
        plt.legend(loc='upper left', fontsize='large')
        plt.grid(color='gray', linestyle=':', linewidth=0.5)
        plt.xlabel("Velocity [m/s]")
        plt.ylabel("Power [kW]")
        plt.show()

    def final_power(self):
        self.P = {
            'HIGE1': self.P_hoge * self.hige_factor,
            'V_climb': 85000, #needs to be implemented
            'HOGE1': self.P_hoge,
            'Climb1': self.min_power_CD, # 9deg climb
            'Cruise1': self.min_power,
            'Descent1': 40000, # steep descent, needs to be implemented
            'HOGE2': self.P_hoge,
            'Loiter': self.min_power,
            'Climb2': 100000, # needs to be implemented
            'Cruise2': self.min_power,
            'Descent2': 40000, # needs to be implemented
            'HOGE3': self.P_hoge,
            'V_Descent': 37500, # needs to be implemented
            'HIGE2': self.P_hoge * self.hige_factor
        }

class SoundAnalysis:
    def __init__(self):
        #Conversions
        self.m_to_f = 1/0.3048 #Conversion factor from meter to feet
        self.N_to_lbs = (1/9.80665)*2.20462

        #Input
        self.rotorsizing = RotorSizing()

        #Observer position relative to rotor
        self.x = 0 #np.linspace(0,100,1000) #Measured in direction of motion of helicopter [ft]
        self.y = 0 #np.linspace(-100,100,1000) #Measured at 90 deg to x in disc plane [ft]
        self.z = 500 #Flyover height [ft]
        self.r = (self.x**2+self.y**2+self.z**2)**0.5

        #Inputs for rotational noise
        self.R = self.rotorsizing.rotor_radius*self.m_to_f #Get rotor radius in [f]
        self.A = np.pi*(self.R**2) #Rotor area [ft^2]
        self.n = self.rotorsizing.omega #Rotor rotational speed [rad/s] 
        self.V = self.rotorsizing.V_max*self.m_to_f #Set to max speed for now [ft/s]
        self.c = self.rotorsizing.speed_of_sound*self.m_to_f #Speed of sound [ft/s]
        self.B = self.rotorsizing.n_blades
        self.T = self.rotorsizing.T_forward_flight*self.N_to_lbs/self.rotorsizing.N_rotors #Thrust [lbs]

        #Inputs for vortex noise
        self.D = 2*self.R #Rotor diameter [ft]
        self.V_07 = 0.7*((self.n*np.pi*self.D)/60) #Linear speed of 0.7 radius section
        self.chord = self.rotorsizing.chord*self.m_to_f #Chord [ft]
        self.A_b = self.B*self.R*self.chord #Total blade area [ft^2]

        #Initialize
        self.rotational_noise()
        self.vortex_noise()

    def rotational_noise(self):
        #Calculate rotational Mach number
        self.M = (0.8*self.n*self.R)/(self.c)

        #Calculate flight Mach number
        self.M_f = self.V/self.c

        #Calculate theta dash
        self.theta_dash = np.arccos(self.x/self.r)

        #Calculate effective rotational Mach number
        self.M_E = self.M/(1-self.M_f*np.cos(self.theta_dash))

        #Calculate theta
        self.theta = np.pi/2 #Set to angle that gives max noise for now (unused)

        #Interpolate from graph of first harmonic
        M_series = np.linspace(0.2,1,9)
        noise_series = np.array([77,82,85,88,91,93,95,96,97])
        interpolator = np.polyfit(M_series, noise_series,3)

        self.rotational_SPL_uncorrected = interpolator[0]*(self.M_E)**3  + interpolator[1]*(self.M_E)**2 + interpolator[2]*(self.M_E) + interpolator [3]

        #Apply SPL correction
        self.rotational_SPL = self.rotational_SPL_uncorrected + 11 + 10*np.log10((self.T/(self.r**2))*(self.T/self.A))

        #Calculate fundamental frequency
        self.f_rotational = (self.n*self.B)/(2*np.pi*(1-self.M_f*np.cos(self.theta)))

    def vortex_noise(self):
        #Calculate SPL for 300 ft distance
        self.vortex_SPL_uncorrected = 10*(2*np.log10(self.V_07)+2*np.log10(self.T)-np.log10(self.A_b)-3.57)

        #Correct for distance (500 ft)
        self.vortex_SPL = self.vortex_SPL_uncorrected - 20*np.log10(self.z/300)

        #Calculate frequency
        self.thickness = 0.12*self.chord #NACA0012 airfoil thickness to chord
        self.f_vortex = (self.V_07*0.28)/self.thickness  #Only valid for small AoA

    def display_parameters_rotor(self):
        print(f"Equivalent M: {self.M_E :.4f}")
        print(f"Uncorrected SPL: {self.rotational_SPL_uncorrected:.4f} dB")
        print(f"Correction factor: {self.rotational_SPL-self.rotational_SPL_uncorrected:.4f} dB")
        print(f"Corrected SPL: {self.rotational_SPL:.4f} dB")
        print(f"Fundamental frequency: {self.f_rotational:.4f} Hz")
        
    def display_paramenters_vortex(self):
        print(f"Uncorrected SPL: {self.vortex_SPL_uncorrected:.4f} dB")
        print(f"Corrected SPL: {self.vortex_SPL:.4f} dB")
        print(f"Vortex frequency: {self.f_vortex:.4f} Hz")


class EnergyAnalysis:
    def __init__(self, V_roc=0.76, rho=1.225, eta_p=0.9, V_climb1=30, V_cruise=40, V_loiter=6, V_descent=35, V_climb2=32, Volts=840, motor_eta=0.85, climb_angle = 9, descent_angle = 5):
        # Assigning the input parameters
        self.Vroc = V_roc
        self.eta_p = eta_p
        self.rho = rho
        self.Volts = Volts
        self.V_climb1 = V_climb1
        self.V_cruise = V_cruise
        self.V_loiter = V_loiter
        self.V_descent = V_descent
        self.V_climb2 = V_climb2
        self.motor_eta = motor_eta  # motor efficiency i think
        self.gamma_1 = climb_angle  # deg
        self.gamma_2 = descent_angle
        
        # Initialize the Power object
        self.power = Power()
        
        # Create a dictionary to hold times, powers, energies, and amps
        self.mission_data = {
            'times': { 'HIGE1': None, 'V_climb': None, 'HOGE1': None, 'Climb1': None, 'Cruise1': None, 'Descent1': None, 'HOGE2': None, 
                      'Loiter': None, 'Climb2': None, 'Cruise2': None, 'Descent2': None, 'HOGE3': None, 'V_Descent': None, 'HIGE2': None, 'total': None },
            'powers': self.power.P,
            'energies': { 'HIGE1': None, 'V_climb': None, 'HOGE1': None, 'Climb1': None, 'Cruise1': None, 'Descent1': None, 'HOGE2': None, 
                          'Loiter': None, 'Climb2': None, 'Cruise2': None, 'Descent2': None, 'HOGE3': None, 'V_Descent': None, 'HIGE2': None, 'total': None },
            'amps': { 'HIGE1': None, 'V_climb': None, 'HOGE1': None, 'Climb1': None, 'Cruise1': None, 'Descent1': None, 'HOGE2': None, 
                      'Loiter': None, 'Climb2': None, 'Cruise2': None, 'Descent2': None, 'HOGE3': None, 'V_Descent': None, 'HIGE2': None, 'max': None }
        }
    
    def calculate_missionphase_time(self):
        # Define given constants
        h1, h2, h3 = 60.0, 300.0, 30.0  # Altitudes (m)
        delta_h2 = h2 - h1
        delta_h3 = h2 - h3
        d_tothalf = 30000  # Distance to half
        
        # Calculate times using loops
        times_dict = self.mission_data['times']
        
        # HIGE 1 time
        times_dict['HIGE1'] = 15
        
        # Vertical Climb time (T2)
        times_dict['V_climb'] = h1 / self.Vroc
        
        # HOGE 1 time
        times_dict['HOGE1'] = 10
        
        # Steady Climb 1 time (T4)
        d_cl1 = delta_h2 / np.tan(np.radians(self.gamma_1))
        times_dict['Climb1'] = d_cl1 / self.V_climb1
        
        # Descent 1 time (T6)
        roc_descent = 7.6 # m/s
        times_dict['Descent1'] = delta_h3 / roc_descent

        # Cruise 1 time (T5)
        d_cr1 = d_tothalf - d_cl1 - (self.V_climb1 * times_dict['Descent1'])
        times_dict['Cruise1'] = d_cr1 / self.V_cruise
        
        # HOGE 2 time
        times_dict['HOGE2'] = 30
        
        # Loiter time
        times_dict['Loiter'] = 600 # Dependent on max endurance
        
        # Steady Climb 2 time (T9)
        d_cl2 = delta_h3 / np.tan(np.radians(self.gamma_1))
        times_dict['Climb2'] = d_cl2 / self.V_climb2
        
        # Descent 2 time (T11)
        d_des2 = delta_h2 / np.tan(np.radians(self.gamma_2))
        times_dict['Descent2'] = d_des2 / self.V_descent

        # Cruise 2 time (T10)
        d_cr2 = d_tothalf - d_cl2 - (self.V_descent * times_dict['Descent2'])
        times_dict['Cruise2'] = d_cr2 / self.V_cruise
        
        # HOGE 3 time
        times_dict['HOGE3'] = 10
        
        # Vertical Descent time (T13)
        Vrod = 0.5
        times_dict['V_Descent'] = h1 / Vrod
        
        # HIGE 2 time
        times_dict['HIGE2'] = 15
        
        # Total time
        # Initialize sum to 0
        total_time = 0
        # Iterate through the 'times' dictionary and sum values excluding 'total'
        for key, value in self.mission_data['times'].items():
            if key != 'total' and value is not None:
                total_time += value
    
        #total_time = sum(times_dict.values()[:-1]) 
        times_dict['total'] = total_time
        
        # Update mission data
        self.mission_data['times'] = times_dict

        return self.mission_data['times']
    
    def calculate_energy_required(self):
        energies_dict = self.mission_data['energies']
        times_dict = self.mission_data['times']
        
        for phase, time in times_dict.items():
            if phase != 'total':
                # Calculate energy (Wh)
                power = self.power.P.get(phase, 0)
                energies_dict[phase] = (time / 3600) * power / self.motor_eta
        
        # Calculate total energy
        # Initialize sum to 0
        total_energy = 0

        # Iterate through the 'energies' dictionary and sum values excluding 'total'
        for key, value in self.mission_data['energies'].items():
            if key != 'total' and value is not None:
                total_energy += value

        energies_dict['total'] = total_energy
        
        # Update mission data
        self.mission_data['energies'] = energies_dict

        return self.mission_data['energies']
    
    def calculate_amps(self):
        amps_dict = self.mission_data['amps']
        powers_dict = self.power.P
        
        for phase, power in powers_dict.items():
            amps_dict[phase] = power / (self.Volts * self.motor_eta)
        
        #amps_dict['max'] = max(amps_dict.values())
        # Initialize max_amps to store the maximum value
        max_amps = None

        # Iterate through the 'amps' dictionary, excluding 'max' and None values
        for key, value in self.mission_data['amps'].items():
            if key != 'max' and value is not None:
                if max_amps is None or value > max_amps:
                    max_amps = value

        amps_dict['max'] = max_amps
        # Update mission data
        self.mission_data['amps'] = amps_dict

        return self.mission_data['amps']
