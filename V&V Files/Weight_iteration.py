import numpy as np 
from math import * 

class weight_iteration:
    def __init__(self, old_battery_w, old_PEMFC_w, old_propulsion_w, battery_energy, fuelcell_power,  omega, chord, blade_n, R,  rotor_n = 1, power_per_cell_avg = 0.181902826055829, hydrogen_weight = 5,  mtow = 718 ):
        #Starting mtow 
        self.mtow = mtow 

        self.omega = omega
        self.chord = chord
        self.blade_n = blade_n
        self.R = R 
        
        self.old_battery_w =  old_battery_w
        self.old_control_devices_w =  0
        self.old_h2_h2tank_w = 73.6
        self.old_PEMFC_w = old_PEMFC_w #ibrahim
        self.old_propulsion_w = old_propulsion_w
        self.old_hydrogen_w = hydrogen_weight
        self.old_control_devices_w = 0 

        self.fuelcell_power = fuelcell_power 
        self.power_per_cell_avg = power_per_cell_avg
        self.battery_energy = battery_energy
        self.hydrogen_mass = hydrogen_weight 
        self.rotor_n = rotor_n
        self.battery_energy = battery_energy

        self.reg_coef_fc = 0.739
        self.fc_offset = 132
        self.reg_coef_battery = 0.00421
        self.battery_offset = 6.95
        self.reg_coef_tank_700 = 13.4
        self.tank_700_offset = 13.3
        self.ref_coef_tank_350 = 9.99
        self.tank_350_offset = 16.6
        
        # variables to be iterated
        self.h2_tank_w = None 
        self.battery_w = None 
        self.PEMFC_w = None 
        self.rotor_w = None 
        self.motor_w = None 

    def calc_new_fuelcell(self):
        self.PEMFC_w = self.fuelcell_power * self.reg_coef_fc + self.fc_offset
        print("PEMFC weight is:", self.PEMFC_w)
        return self.PEMFC_w
    
    def calc_new_battery(self):
        self.battery_w = self.battery_energy * self.reg_coef_battery + self.battery_offset
        print("battery weight is", self.battery_w)
        return self.battery_w
    
    #weight difference is minor in comparison to the 
    def calc_new_tank(self):
        self.h2_tank_w = self.reg_coef_tank_700 * self.hydrogen_mass + self.tank_700_offset 
        self.h2_tank_w
        return self.h2_tank_w

    def calc_new_rotor_w(self):
        self.rotor_w = 0.026 * (self.blade_n ** 0.66) * self.chord * (self.R ** 1.3) * (self.omega * self.R)**0.67 
        print("rotor weight is", self.rotor_w)
        return self.rotor_w 
    
    def calc_new_motor_w(self):
        self.motor_w = float(input("what is the new motor weight you found off the shelf?"))
        return self.motor_w

    def update_new_mtow(self): 
        self.calc_new_battery()
        self.calc_new_fuelcell()
        self.calc_new_rotor_w()
        self.calc_new_tank()
        self.calc_new_motor_w()
        new_mtow = self.mtow - self.old_battery_w - self.old_control_devices_w - self.old_h2_h2tank_w - self.old_PEMFC_w - self.old_propulsion_w + self.PEMFC_w + self.battery_w + self.h2_tank_w + self.rotor_w + self.motor_w
        #if new_mtow > self.mtow:
        self.mtow = new_mtow
        return self.mtow 
    
    


coax_w = weight_iteration(old_battery_w=98.187, old_PEMFC_w=196.39, old_propulsion_w= 51.282, battery_energy= 15710, fuelcell_power=87.139, omega= 504.01, chord=0.198, blade_n=2, R=3.89)
print ("coax_iterated_weight is :", coax_w.update_new_mtow())

quad_W = weight_iteration(old_battery_w=82.872, old_PEMFC_w=191.41, old_propulsion_w= 24.13, battery_energy= 13259 , fuelcell_power= 80.405, omega=1279.275, chord = 0.15, blade_n=4, R = 1.29)
print("quad_iterated_weight is :",quad_W.update_new_mtow())

compound_w = weight_iteration(old_battery_w=379.591, old_PEMFC_w=147.28, old_propulsion_w= 131.086, battery_energy= 58633, fuelcell_power=20.689 ,omega=3400, chord = 0.254,blade_n=5, R = 0.98 )
print("compound_iterated_weight is :", compound_w.update_new_mtow())
3