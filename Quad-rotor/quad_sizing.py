import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

class RotorAnalysis:
    def __init__(self, mtow_kg=718.89, n_rotors=4, dl_imperial=7, bank_angle=30, v_max_kmh=150):
        # Constants
        self.MTOW_KG = mtow_kg  # Maximum Takeoff Weight in kg
        self.N_rotors = n_rotors  # Number of rotors
        self.DL_IMPERIAL = dl_imperial  # Disc loading in lb/ft²
        self.bank_angle = bank_angle  # Bank angle during turn [deg]
        self.V_MAX_KMH = v_max_kmh  # Maximum forward speed in km/h

        # Physical constants
        self.POUNDS_PER_KG = 2.205  # Conversion factor from kg to pounds
        self.FT_TO_M = 3.281  # Conversion factor from ft to m
        self.SPEED_OF_SOUND = 343  # Speed of sound in m/s
        self.rho = 1.225  # Air density [kg/m³]
        self.g = 9.80665  # Gravitational acceleration [m/s²]

        # Derived values
        self.D_v = 0.05 * self.MTOW_KG  # Drag penalty [kg]
        self.k_dl = (1 + (self.D_v / self.MTOW_KG))  # Drag load factor
        self.V_ne = (self.V_MAX_KMH / 3.6) * 1.1  # Never exceed speed [m/s]

        self.V_gust = 30 / self.FT_TO_M
        self.C_lalpha = 5.73  # 1/rad; NACA0012

        #User inputs 
        self.c_t_o_fl = 1 # random, gets updated
        self.c_t_o_turn = 1 # random, gets updated
        self.c_t_o_turb = 1 # random, gets updated

        # Updated variables
        self.radius_m = None
        self.rotor_diameter_m = None
        self.tip_speed = None
        self.mach_num = None
        self.omega = None
        self.adv_ratio_fl = None
        self.n_z = None
        self.o_fl = None
        self.o_turn = None
        self.o_turb = None
        self.solidity = None
        self.aspect_ratios = None
        self.N_bl = None
        self.chord = None
        self.max_c_t = None


    def calculate_radius(self):
        mtow_pounds = self.MTOW_KG * self.POUNDS_PER_KG
        radius_ft = np.sqrt((mtow_pounds / self.N_rotors) / (self.DL_IMPERIAL * np.pi))
        radius_m = radius_ft / self.FT_TO_M
        self.radius_m = radius_m
        return self.radius_m

    def calculate_tip_speed(self):
        tip_speed = 603 / 3.6  # Tip speed in m/s (fixed)
        self.tip_speed = tip_speed
        return self.tip_speed

    def calculate_mach_number(self):
        self.calculate_tip_speed()
        mach_num = ((self.V_MAX_KMH / 3.6) + self.tip_speed) / self.SPEED_OF_SOUND
        self.mach_num = mach_num
        return self.mach_num

    def perform_analysis(self):
        self.rotor_radius_m = self.calculate_radius()
        self.rotor_diameter_m = 2 * self.rotor_radius_m

        # Tip speed
        self.tip_speed = self.calculate_tip_speed()
        self.mach_num = self.calculate_mach_number()

        # Rotational speed
        self.omega = self.tip_speed / self.rotor_radius_m  # Rotational speed in rad/s

        # Advance ratio
        self.adv_ratio_fl = self.V_ne / (self.omega * self.rotor_radius_m)

        # Thrust and solidity in forward flight
        T_fl = self.k_dl * self.MTOW_KG * self.g
        c_t_fl = T_fl / (self.N_rotors * (self.rho * np.pi * self.rotor_radius_m**2 * (self.omega * self.rotor_radius_m)**2))
        self.o_fl = c_t_fl / self.c_t_o_fl

        # Solidity during turn
        self.n_z = 1 / np.cos(self.bank_angle * (np.pi / 180))  # Load factor
        T_turn = self.n_z * self.k_dl * self.MTOW_KG * self.g
        c_t_turn = T_turn / (self.N_rotors * (self.rho * np.pi * self.rotor_radius_m**2 * (self.omega * self.rotor_radius_m)**2))
        self.o_turn = c_t_turn / self.c_t_o_turn

        # Solidity during turbulence
        delta_n = (0.25 * self.C_lalpha * (self.V_gust / (self.omega * self.rotor_radius_m))) / self.c_t_o_turb
        n_z_turb = 2 + delta_n  # 2g pull up (FAA requirement)
        T_turb = n_z_turb * self.k_dl * self.MTOW_KG * self.g
        c_t_turb = T_turb / (self.N_rotors * (self.rho * np.pi * self.rotor_radius_m ** 2 * (self.omega * self.rotor_radius_m) ** 2))
        self.o_turb = c_t_turb / self.c_t_o_turb

        # Maximum solidity
        self.solidity = max(self.o_fl, self.o_turn, self.o_turb)
        self.max_c_t = max(c_t_fl, c_t_turn, c_t_turb)

        # Aspect ratio calculation for different numbers of blades
        N_blades = range(2, 7)
        self.aspect_ratios = []

        for self.N_bl in N_blades:
            chord = (self.solidity * np.pi * self.rotor_radius_m) / self.N_bl
            self.AR_blades = self.rotor_radius_m**2 / (self.rotor_diameter_m * chord)
            self.aspect_ratios.append((self.N_bl, chord, self.AR_blades))

        # Return all results including advance ratio
        return self.rotor_radius_m, self.rotor_diameter_m, self.tip_speed, self.mach_num, self.omega, self.adv_ratio_fl, self.o_fl, self.o_turn, self.o_turb, self.solidity, self.aspect_ratios, self.n_z

    def plot_aspect_ratios(self):
        N_blades = [x[0] for x in self.aspect_ratios]
        aspect_ratios_values = [x[2] for x in self.aspect_ratios]

        plt.figure(figsize=(8, 5))
        plt.plot(N_blades, aspect_ratios_values, marker='o', label="Aspect Ratio")
        plt.title("Aspect Ratio vs. Number of Blades")
        plt.xlabel("Number of Blades")
        plt.ylabel("Aspect Ratio")
        plt.grid(True)
        plt.legend()
        #plt.show()

    def calculate_chord(self):
        """Calculate the blade chord length given the number of blades and solidity."""
        chord = (self.solidity * np.pi * self.rotor_radius_m) / self.N_bl
        self.chord = chord
        return self.chord

    def interact(self):
        # Perform analysis with default coefficients to get advance ratio first
        results = self.perform_analysis()

        # Print advance ratio and other results
        print(f"Advance Ratio in forward flight = {self.adv_ratio_fl:.2f}")
        #print(f"Rotor Radius [m]: {self.rotor_radius_m:.2f}")
        #print(f"Rotor Diameter [m]: {self.rotor_diameter_m:.2f}")
        #print(f"Tip Speed [m/s]: {self.tip_speed:.2f}")
        #print(f"Mach Number: {self.mach_num:.2f}")
        #print(f"Rotational Speed (Omega) [rad/s]: {self.omega:.2f}")
        #print(f"RPMs : {self.omega*(60/(2*np.pi)): .2f}")

        # Now ask for user input for thrust coefficients
        #print("Please enter the thrust coefficients:")
        self.c_t_o_fl = float(input("Enter thrust coefficient for forward flight (c_t_o_fl): "))
        self.c_t_o_turn = float(input("Enter thrust coefficient for turning flight (c_t_o_turn): "))
        self.c_t_o_turb = float(input("Enter thrust coefficient for turbulence (c_t_o_turb): "))

        # Perform analysis with user-provided coefficients
        results_with_input = self.perform_analysis()

        # Print new results with user inputs
        #print(f"Solidity in Forward Flight: {self.o_fl:.4f}")
        #print(f"Solidity during Turn: {self.o_turn:.4f}")
        #print(f"Solidity during Turbulence: {self.o_turb:.4f}")
        #print(f"Solidity used: {self.solidity:.4f}")

        # Plot results
        #self.plot_aspect_ratios()

        # Ask user for the number of blades they want and calculate the corresponding chord
        selected_blades = int(input("Enter the number of blades you want (between 2 and 6): "))
        if selected_blades not in range(2, 7):
            raise ValueError("Invalid number of blades selected. Please choose between 2 and 6 blades.")
        else:
            for self.N_bl, self.chord, self.AR_blades in self.aspect_ratios:
                if self.N_bl == selected_blades:
                    self.AR_blades = self.AR_blades
                    self.chord = self.chord

            # Print the results
            #print(f"Selected number of blades: {selected_blades}")
            #print(f"Corresponding blade chord length: {self.chord:.4f} meters")
            #print(f"Corresponding aspect ratio: {self.AR_blades:.4f}")


class PerformanceAnalysis:
    def __init__(self, rotor_analysis: RotorAnalysis, vertical_climb_speed: float = 0.76):
        self.rotor_analysis = rotor_analysis
        self.V_vertical_climb = vertical_climb_speed

        # Use values from RotorAnalysis
        self.mtow = rotor_analysis.MTOW_KG
        self.N_rotors = rotor_analysis.N_rotors
        self.rho = rotor_analysis.rho
        self.g = rotor_analysis.g
        self.Cl_alpha = rotor_analysis.C_lalpha
        self.mtowN = self.mtow * self.g

        # Perform rotor sizing
        self.results = rotor_analysis.perform_analysis()
        self.rotor_radius, self.rotor_diameter, self.v_tip, self.mach_number, self.omega, self.adv_ratio_fl, self.o_fl, self.o_turn, self.o_turb, self.solidity, self.aspect_ratios, self.n_z = self.results

        # Derived parameters
        self.rotor_area = np.pi * (self.rotor_radius ** 2)  # Rotor swept area (m^2)
        self.weight_per_rotor = (self.mtow / self.N_rotors) * self.g  # Weight supported by each rotor (N)
        self.FM = 0.78  # figure of merit (assumed)

        # Induced power calculation constants
        self.V_point = 13  # cruise speed, assumption for now
        self.T = rotor_analysis.k_dl * self.mtowN
        self.k = 1.15  # inflow correction, assumption: should be between 1.1 and 1.2

    def calculate_hover_induced_velocity(self):
        """Calculate hover induced velocity (m/s)."""
        return np.sqrt(self.weight_per_rotor / (2 * self.rho * self.rotor_area))

    def calculate_climb_induced_velocity(self):
        """Calculate climb induced velocity (m/s)."""
        return (-self.V_vertical_climb / 2) + np.sqrt((self.V_vertical_climb / 2)**2 + (self.weight_per_rotor / (2 * self.rho * self.rotor_area)))

    def calculate_disk_area(self):
        """Calculate total disk area of the quadcopter"""
        S_disk = self.rotor_area * 4
        return S_disk

    def calculate_power_req_hover(self):
        """Calculate the required power for hover."""
        P_req_hover = ((self.mtow * self.g) ** (3 / 2)) / (np.sqrt(2 * self.rho * self.calculate_disk_area()) * self.FM)
        return P_req_hover
    
    def calculate_power_loading_hover(self):
        P_load_hover = (self.mtow) / (self.calculate_power_req_hover()/1000)
        return P_load_hover

    def calculate_C_T(self):
        """Calculate the thrust coefficient."""
        C_T = self.T / (self.rho * self.calculate_disk_area() * (self.omega * self.rotor_radius) ** 2)
        return C_T

    def calculate_average_C_l(self):
        """Calculate average coefficient of lift."""
        average_C_l = 6.6 * self.calculate_C_T() / self.solidity
        return average_C_l

    def calculate_alpha_m(self):
        """Calculate angle of attack for maximum lift."""
        alpha_m = self.calculate_average_C_l() / self.Cl_alpha
        return alpha_m

    def calculate_average_C_D_p(self):
        """Calculate the average profile drag coefficient."""
        average_C_D_p_1 = 0.0087 - 0.0216 * self.calculate_alpha_m() + 0.4 * self.calculate_alpha_m() ** 2  # bailey
        average_C_D_p_2 = 0.011 + 0.4 * (self.calculate_alpha_m() ** 2)  # marinescu
        average_C_D_p_3 = 0.009 + 0.73 * (self.calculate_alpha_m() ** 2)  # talbot
        # pick the highest
        average_C_D_p = max(average_C_D_p_1, average_C_D_p_2, average_C_D_p_3)
        return average_C_D_p

    def calculate_P_p_hov(self):
        """Calculate the profile drag power for hover."""
        P_p_hov = self.solidity * self.calculate_average_C_D_p() / 8 * self.rho * (self.omega * self.rotor_radius) ** 3 * self.calculate_disk_area()
        return P_p_hov

    def calculate_P_p(self):
        """Calculate the profile drag power for forward flight."""
        P_p = self.calculate_P_p_hov() * (1 + 4.65 * self.adv_ratio_fl ** 2)
        return P_p

    def calculate_v_i_bar(self):
        """Equation to solve for induced velocity ratio."""
        v_i_hover = self.calculate_hover_induced_velocity()
        def equation(v_i_bar):
            return v_i_bar ** 4 + (self.V_point / v_i_hover) ** 2 * v_i_bar ** 2 - 1
        initial_guess = 1.0
        solution = fsolve(equation, initial_guess)
        return solution[0]

    def calculate_P_i(self):
        """Calculate induced power."""
        # a little uncertainty here from me - Jan
        # we use the induced velocity of one rotor and then multiply it by the total thrust getting power induced
        # this kinda confuses me but if i divide thurst by 4 (each rotor induced power) and then multiply power by 4 at the end i get the same result
        v_i_fl = self.calculate_v_i_bar() * self.calculate_hover_induced_velocity()
        P_i = self.k * (self.T/4) * v_i_fl * 4
        return P_i

    def calculate_P_par(self):
        """Calculate parasite power."""
        A_eq = 0.418  # m^2, assuming clean configuration
        P_par = 0.5 * A_eq * self.rho * (self.V_point ** 3)
        return P_par

    def calculate_P_req_fl(self):
        """Calculate the power required for forward flight."""
        P_req_fl = (self.calculate_P_p() + self.calculate_P_i() + self.calculate_P_par()) * 1.045  # accounting for 3% to 6% losses
        return P_req_fl

    def display_results(self):
        """Display power requirements and induced velocities."""
        v_i_hover = self.calculate_hover_induced_velocity()
        v_i_climb = self.calculate_climb_induced_velocity()
        disk_area = self.calculate_disk_area()
        power_requirements = self.calculate_power_req_hover()
        power_loading_hover = self.calculate_power_loading_hover()
        profile_power_fl = self.calculate_P_p()
        power_required_forward = self.calculate_P_req_fl()

        print(f'Disk Area: {disk_area:.2f} m²')
        print(f'Power required for hover: {power_requirements / 1000:.2f} kW')
        print(f'Power loading for hover: {power_loading_hover:.2f} kg/kW')
        print(f'Power required for forward flight (not optimal speed): {power_required_forward / 1000:.2f} kW')

    def calculate_power_curves(self):
        """Calculate power requirements for a range of forward speeds and display results."""
        # Define the range of forward speeds
        V_points = np.linspace(5, 60, 1000)  # Forward speeds from 5 to 60 m/s

        # Arrays to store solutions
        P_is = []  # Induced power
        P_pars = []  # Parasite power
        P_reqs = []  # Total required power

        # Calculate power values for each forward speed
        for V_point in V_points:
            self.V_point = V_point  # Update forward speed
            P_is.append(self.calculate_P_i() / 1000)  # Induced power in kW
            P_pars.append(self.calculate_P_par() / 1000)  # Parasite power in kW
            P_reqs.append(self.calculate_P_req_fl() / 1000)  # Total required power in kW

        # Find the minimum power and corresponding speed
        min_power = min(P_reqs)
        min_speed = V_points[P_reqs.index(min_power)]

        print(f'The minimum power required in forward flight is: {min_power:.2f} kW '
              f'and the corresponding optimal speed is: {min_speed * 3.6:.2f} km/h')

        # Plot the power curves
        plt.figure(figsize=(12, 8))
        plt.plot(V_points * 3.6, P_is, label='Induced Power', color='green', linestyle='--', linewidth=2)
        plt.plot(V_points * 3.6, P_pars, label='Parasite Power', color='red', linestyle='-.', linewidth=2)
        plt.plot(V_points * 3.6, P_reqs, label='Total Required Power', color='blue', linewidth=2)

        # Add a point for the minimum power and corresponding speed
        plt.scatter(min_speed * 3.6, min_power, color='purple', zorder=5,
                    label=f'Minimum Power: {min_power:.1f} kW at {min_speed * 3.6:.1f} km/h')

        # Add vertical and horizontal lines to the point
        plt.axvline(x=min_speed * 3.6, color='purple', linestyle=':', linewidth=1)
        plt.axhline(y=min_power, color='purple', linestyle=':', linewidth=1)
        plt.axhline(y=self.calculate_power_req_hover() / 1000, label='Hover Required Power', color='orange', linewidth=2)

        # Add labels, title, legend, and grid
        plt.title('Power Requirements for Forward Flight', fontsize=16)
        plt.xlabel('Forward Speed, $V_{point}$ (km/h)', fontsize=14)
        plt.ylabel('Power (kW)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.show()


"""
class PowerAnalysis:
    def __init__(self, V_roc = 0.76, eta_p = 0.9):
        #values from rotor sizing
        self.N_rotors = RotorAnalysis.N_rotors     # get the amount of rotors    
        #values from performance
        self.P_req_hover = PerformanceAnalysis.calculate_power_req_hover()     # get the required power for hover
        self.mtowN = PerformanceAnalysis.mtowN
       
        self.Vroc = V_roc       #Vertcial rate of climb
        self.eta_p = eta_p      #propulsion efficiency
        
        #add power required for take off, climb, cruise, descent, loiter, land
        self.t4 = None
        #NEED TO BE ADDED - import the V_climb1 FROM PERFORMANCE!! - Ask Yunjae!
        #V_climb1 = 36.0       #[m/s] test value


    
    def calculate_missionphase_time(self):
        # HIGE 1
        t1 = 15                #[s] - Given
        # Vertical Climb
        delta_h1 = 60.0     #[m] - Given
        t2 = delta_h1 / self.Vroc
        # HOGE 1  
        t3 = 10             #[s] - Given
        # Steady Climb 1
        h2 = 300        #[m] - Given
        delta_h2 = h2-delta_h1      #[m] height between Hoge1 and cruise 1
        gamma_1 = 9         #[deg] flight path angle
        d_cl1 = delta_h2 / np.tan(gamma_1*np.pi/180)      #Horizontal distance after converting deg to rad
        t4 = d_cl1 / V_climb1        #time of climb 1 with horizontal climb speed V_climb1
        # Descent 1
        h3 = 30         #[m]       Hoge 2 altitude 
        roc_descent = 7.6      #[m/s]   -    given
        delta_h3 = h2 - h3      #[m] height between cruise 1 and Hoge2
        t6 = delta_h3 / roc_descent     #[s]
        d_ds1 = V_climb1 * t6    # horizontal distance covered during descent
        # Cruise 1
        d_tothalf = 30000       #[m] - Given distance from wright memorial to alligators
        d_cr1 = d_tothalf - d_cl1 - d_ds1       #[m] Horizontal cruise1 distance
        t5 = d_cr1 / V_cruise       #[s] time needed for cruise1, depending on the cruise speed that has to be implimented!
        # HOGE 2
        t7 = 30     #[s]
        # Loitering
        t8 = t8     #DEPENDS ON MAX ENDURENCE!!-> to be designed for
        d_loiter = V_loiter * t8
        # Steady Climb 2
        delta_h4 = h2 - h3
        gamma_2 = 9     #[deg] flight path angle
        d_cl2 = delta_h4 / np.tan(gamma_2*np.pi/180)        #[m] horizontal distance of climb 2
        t9 = d_cl2 / V_climb2       #[s] 
        # Descent 2
        delta_h5 = h2-delta_h1
        gamma_3 = 5     #[deg] flight path angle down
        d_ds2 = delta_h5 / np.tan(gamma_3*np.pi/180)        # horizontal distance covered during descent 2
        t11 = d_ds2 / V_descent     #[s] time to cover horizontal distance according to horizontal descent speed
        # Cruise 2
        d_cr2 = d_tothalf - d_cl2 - d_ds2       #[m] Horizontal cruise2 distance
        t10 = d_cr2 / V_cruise       #[s] time needed for cruise2, depending on the cruise speed that has to be implimented!
        # HOGE 3
        t12 = 10        #[s]
        # Vertical Descent
        Vrod = 0.5      #[m/s]  Vertical rate of descent
        t13 = delta_h1 / Vrod       #[s] Vertical descent time
        # HIGE 2
        t14 = 15        #[s]

        t_total = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13 + t14
        d_total = 2 * d_tothalf + d_loiter

        return t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t_total, d_total
    
    def calculate_power_required(self):

        # power required for HOGE
        p_hover = self.P_req_hover   #import P_req_hover from PerformanceAnalysis

        # power required for vertical climb
        P_vclimb = self.Vroc * self.mtowN     #import Vroc from MTOW in Newtons from PerformanceAnalysis
        P_vclimbeff = P_vclimb / self.eta_p        # eta_p is the propulsion efficiency
"""


if __name__ == "__main__":
    # Initialize the analysis with default values
    rotor_analysis = RotorAnalysis()

    # Start the interaction method that handles all the user input and analysis
    rotor_analysis.interact()

    # Initialize PerformanceAnalysis with the rotor analysis instance
    performance_analysis = PerformanceAnalysis(rotor_analysis)

    # Display results for hover and basic performance
    #performance_analysis.display_results()

    # Calculate and display power curves for forward flight
    #performance_analysis.calculate_power_curves()