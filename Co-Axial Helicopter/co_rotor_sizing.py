import math
import numpy as np
import matplotlib.pyplot as plt

class RotorSizing:

    def __init__(self, MTOW=718.89, n_blades=6):

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

        # initialize the parameters 
        self.update_parameters()

    def update_parameters(self):
        # derived parameters 
        self.disc_loading = 2.28 * ((self.MTOW)**(1/3) - 2.34)                      # kg/m^2
        self.rotor_radius = math.sqrt(self.MTOW / (math.pi * self.disc_loading))    # m 
        self.rotor_diameter = 2 * self.rotor_radius                                 # m
        self.tip_speed = 140 * (self.rotor_diameter)**0.171                         # m / s
        self.RPM = 2673 / ((self.rotor_diameter)**0.892)                            # rpm
        self.omega = self.RPM_to_rad * self.RPM                                     # rad / s

        # flight performance 
        self.max_forward_velocity = (self.max_tip_mach * self.speed_of_sound) - self.tip_speed  # m / s
        self.never_exceed_velocity = 1.1 * self.max_forward_velocity                            # m / s
        self.mu_ne = self.never_exceed_velocity / self.tip_speed                                # --

        # forward flight solidity 
        self.fuselage_download = 0.04 * self.MTOW * self.g                 # N
        self.k_dl = 1 + (self.fuselage_download / (self.MTOW * self.g))    # --
        self.T_forward_flight = self.k_dl * self.MTOW * self.g             # N
        self.C_T_forward_flight = self.T_forward_flight / (self.rho * math.pi * self.rotor_radius**2 * (self.rotor_radius * self.omega)**2) # --
        self.solidity_forward_flight = self.C_T_forward_flight / 0.08                     # --

        # turn solidity 
        self.roll_angle = 30                                        # deg
        self.n_z_turn = 1 / np.cos(self.roll_angle * self.deg_to_rad)    # --
        self.T_turn = self.n_z_turn * self.k_dl * self.MTOW * self.g     # N 
        self.C_T_turn = self.T_turn / (self.rho * math.pi * (self.rotor_radius**2) * (self.rotor_radius * self.omega)**2)   # -- 
        self.solidity_turn = self.C_T_turn / 0.12

        # turbulent solidity
        self.lift_slope = 5.73      # 1 / rad (from NACA0012)
        self.gust_velocity = 30 * self.ft_to_m
        self.delta_n_turbulent = 0.25 * self.lift_slope * (self.gust_velocity/(self.rotor_radius * self.omega))
        self.n_z_turbulent = self.delta_n_turbulent + 2
        self.T_turbulence = self.n_z_turbulent * self.k_dl * self.MTOW * self.g
        self.C_T_turbulence = self.T_turbulence / (self.rho * math.pi * (self.rotor_radius**2) * (self.rotor_radius * self.omega)**2)
        self.solidity_turbulent = self.C_T_turbulence / 0.17

        # chord and aspect ratio 
        self.maximum_solidity = np.max([self.solidity_forward_flight, self.solidity_turn, self.solidity_turbulent])
        self.chord = (self.maximum_solidity * math.pi * self.rotor_radius) / (2 * self.n_blades)
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



if __name__ == '__main__':
    rotor = RotorSizing()
    print('----------------------------------------')
    print(f'Number of Blades = {rotor.n_blades} and MTOW = {rotor.MTOW}kg')
    print('----------------------------------------')
    rotor.display_parameters()
    ### visualizations ###
    rotor.visual_blade_vs_aspect_ratio()