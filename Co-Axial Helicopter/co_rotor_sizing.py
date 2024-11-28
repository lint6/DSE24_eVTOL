import math
import numpy as np
import matplotlib as plt

class RotorSizing:

    def __init__(self):

        # ambient constants 
        self.rho_ISA = 1.225        #kg/m^3
        self.g = 9.80665        #m/s^2 
        self.temperature_ISA = 15         #celsius
        self.speed_of_sound_ISA = math.sqrt(1.4 * 287 * (self.celsius_to_kelvin + self.temperature_ISA))



        # physical constants 
        self.MTOW = 718.89          #kg 
        self.n_blades = 6           # -
        self.max_tip_mach = 0.85    # -

        # conversion values
        self.kg_to_N = 9.80665              #multiply
        self.celcius_to_kelvin = 273.15     #add
        self.RPM_to_rad = (math.pi / 30)    #multiply
        self.deg_to_rad = math.pi / 180     #multiply


        # initial calculated values (from statistic relationships)
        self.disc_loading = 2.28 * ((self.MTOW)**(1/3) - 2.34)                      #kg/m^2
        self.rotor_radius = math.sqrt(self.MTOW / (math.pi * self.disc_loading))    #m
        self.rotor_diameter = 2 * self.rotor_radius                                 #m
        self.tip_speed = 140 * (self.rotor_diameter)**0.171                         #m/s
        self.RPM = 2673 / ((self.rotor_diameter)**0.892) 
        self.omega = self.RPM_to_rad * self.RPM                        #m                            
        self.chord = 0.0108 * ((self.MTOW**0.539) / (self.n_blades**0.714))         #m
        self.solidity = (self.n_blades * self.chord) / (math.pi * self.radius)      # -
        self.max_forward_velocity = (self.max_tip_mach * self.speed_of_sound_ISA) - self.tip_speed  #m/s
        self.never_exceed_velocity = 1.1 * self.max_forward_velocity 

        # blade solidity in level-flight 
        self.mu_ne = self.never_exceed_velocity / self.tip_speed
        self.fuselage_download = 0.04 * self.MTOW * self.kg_to_N #N (4% from MTOW, statistics)
        self.k_dl = 1 + (self.fuselage_download / (self.MTOW * self.kg_to_N))
        self.T_forward_flight = self.k_dl * self.MTOW * self.kg_to_N
        self.C_T = self.T_forward_flight / (self.rho_ISA * math.pi * (self.rotor_radius)**2 * (self.rotor_radius * self.omega)**2)
        self.max_ct_to_solidity = 0.08 # from graph
        self.solidity_lf = self.C_T / self.max_ct_to_solidity

        # blade solidity at turn 
        self.roll_angle = 30 #deg (from regulations)
        self.n_z = 1 / np.cos(self.roll_angle * self.deg_to_rad)
        self.T_turn = self.n_z * self.k_dl * self.MTOW * self.kg_to_N
        self.C_T_turn = self.T_turn / (self.rho_ISA * math.pi * (self.rotor_radius)**2 * (self.rotor_radius * self.omega)**2)
        self.max_ct_to_solidity = 0.12 # from graph 
        self.solidity_turn = self.C_T_turn  / self.max_ct_to_solidity
        








