# Made by Ibrahim and Floris

# This tool outputs an initial estimate of where main components should be placed within specified bounds in the aircraft to achieve the desired CG and that can be used to get a fuselage layout. 

# For Quad, inputs are: Mass of fuselage, H2tank, FC, Batteries
# For Coaxial, inputs are: Mass of fuselage, H2tank, FC, Batteries
# For Compound, inputs are: Masses of wing(s) + rotor. Mass of fuselage, H2tank, FC, Batteries

# Please check if the bounds make sense for your given aircraft type!!!

import numpy as np
import matplotlib.pyplot as plt

class MassBalance:

    def __init__(self, aircraft_type):
        self.aircraft_type = aircraft_type
        self.MTOW = 718.89
        self.m_payload = 185
        self.m_fuselage = 100 #ESTIMATE OR INPUT NEW VALUE
        self.m_batteries = 82.872
        self.m_fuelcell = 191.41
        self.m_h2tank = 73.6
        if self.aircraft_type == 1:
            # Define Ideal Total CG for Quadrotor
            self.cg_ideal_ac = 50
            # Define ideal min/max CG positions for each component
            self.payload_cg_min = 10 
            self.payload_cg_max = 15

            self.fuselage_cg_min = 45
            self.fuselage_cg_max = 55

            self.batteries_cg_min = (self.payload_cg_min + 10)
            self.batteries_cg_max = 80

            self.fuelcell_cg_min = (self.payload_cg_min +10)
            self.fuelcell_cg_max = 80
            
            self.h2tank_cg_min = (self.payload_cg_min +10)
            self.h2tank_cg_max = 80
            # Define rotor set parameters
            self.num_rotor_sets = 1
             # Define rotor radius
            self.rotor_radius = 0 #2  # in m
            self.rotor_m = 5 *2 #in kg
            self.rotor_1_cg_min = self.cg_ideal_ac #(self.rotor_radius * 10) 
            self.rotor_1_cg_max = self.cg_ideal_ac #(self.rotor_1_cg_min + 10)
            # self.rotor_2_cg_min = self.cg_ideal_ac-1 #(100- self.rotor_radius -10)
            # self.rotor_2_cg_max = self.cg_ideal_ac #100 - (self.rotor_radius *10)
           
            # Wing parameters, set to 0 for this type of aircraft
            self.wing_1_cg_min = 0
            self.wing_1_cg_max = 0
            self.wing_2_cg_max = 0
            self.wing_2_cg_min = 0
            self.wing_1_m = 0
            self.wing_2_m = 0
            
            self.rotor_cg_positions = [self.rotor_1_cg_min]
        elif self.aircraft_type == 2:
            # Define parameters for Coaxial Helicopter
            self.cg_ideal_ac = 39 
            # Main Components CG bounds
            self.payload_cg_min = 10
            self.payload_cg_max = 15

            self.fuselage_cg_min = 30
            self.fuselage_cg_max = 45

            self.batteries_cg_min = 31
            self.batteries_cg_max = 56

            self.fuelcell_cg_min = 31
            self.fuelcell_cg_max = 56

            self.h2tank_cg_min = 31
            self.h2tank_cg_max = 56
            # Define rotor set parameters
            self.num_rotor_sets = 1
            self.rotor_m = 36 #Based on 3.9 m radius from this site https://aircommand.com/pages/rotor-blade-selection-and-planning#
            self.rotor_1_cg_min = self.cg_ideal_ac 
            self.rotor_1_cg_max = self.cg_ideal_ac 
            self.rotor_cg_positions = [self.cg_ideal_ac]
            # Wing parameters, set to 0 for this type of aircraft
            self.wing_1_cg_min = 0
            self.wing_1_cg_max = 0
            self.wing_2_cg_max = 0
            self.wing_2_cg_min = 0
            self.wing_1_m = 0
            self.wing_2_m = 0
            pass
        elif self.aircraft_type == 3:
            # Define parameters for Compound Aircraft
            self.cg_ideal_ac = 401 #ASK COMPOUND GROUP, lintong said it should be between both ducts but that messes up the code
            self.payload_cg_min = 10
            self.payload_cg_max = 15

            self.fuselage_cg_min = 45
            self.fuselage_cg_max = self.cg_ideal_ac

            self.batteries_cg_min = 30
            self.batteries_cg_max = 70

            self.fuelcell_cg_min = 50
            self.fuelcell_cg_max = 90

            self.h2tank_cg_min = 35
            self.h2tank_cg_max = 70
            # Wing parameters
            self.wing_1_cg_min = (self.payload_cg_min +10)
            self.wing_1_cg_max = (self.payload_cg_max +10)
            self.wing_2_cg_min = (75)
            self.wing_2_cg_max = (95)
            self.wing_1_m = 50 #TO BE CHANGED
            self.wing_2_m = 50 #TO BE CHANGED
             # Define rotor set parameters
            self.num_rotor_sets = 2
            self.rotor_m = 10 #ASK COMPOUND GROUP
            self.rotor_1_cg_min = (self.wing_1_cg_min) #ASK COMPOUND GROUP
            self.rotor_1_cg_max = (self.wing_1_cg_max) #ASK COMPOUND GROUP
            self.rotor_2_cg_min = (self.wing_2_cg_min)
            self.rotor_2_cg_max = (self.wing_2_cg_max)
            self.rotor_cg_positions = [self.rotor_1_cg_min, self.rotor_2_cg_max]
            pass
        else:
            raise ValueError("Invalid aircraft type. Choose 1 for Quadrotor, 2 for Coaxial Helicopter, or 3 for Compound Aircraft.")

        self.cg_payload = np.random.uniform(self.payload_cg_min, self.payload_cg_max)
        self.cg_fuselage = np.random.uniform(self.fuselage_cg_min, self.fuselage_cg_max)
        self.cg_batteries = np.random.uniform(self.batteries_cg_min, self.batteries_cg_max)
        self.cg_fuel_cell = np.random.uniform(self.fuelcell_cg_min, self.fuelcell_cg_max)
        self.cg_h2tank = np.random.uniform(self.h2tank_cg_min, self.h2tank_cg_max)
        self.cg_wing_1 = np.random.uniform(self.wing_1_cg_min, self.wing_1_cg_max)
        self.cg_wing_2 = np.random.uniform(self.wing_2_cg_min, self.wing_2_cg_max)
 
        self.balancing()

    def get_rotorset_m_and_cgs(self):
        m_rotorset = [self.rotor_m] * self.num_rotor_sets * 2
        cgs = []
        for cg in self.rotor_cg_positions:
            cgs.extend([cg] * self.num_rotor_sets)
        return m_rotorset, cgs

    def balancing(self):
        m_rotorset, rotor_cgs = self.get_rotorset_m_and_cgs()
        self.total_m = np.sum([self.m_payload, self.m_fuselage,
                               self.m_batteries, self.m_fuelcell,
                               self.m_h2tank, self.wing_1_m, self.wing_2_m] + m_rotorset)
        self.total_m_weighted = np.sum([self.m_payload * self.cg_payload,
                                        self.m_fuselage * self.cg_fuselage,
                                        self.m_batteries * self.cg_batteries,
                                        self.m_fuelcell * self.cg_fuel_cell,
                                        self.m_h2tank * self.cg_h2tank, 
                                        self.wing_1_m * self.cg_wing_1,
                                        self.wing_2_m * self.cg_wing_2] +  [m * cg for m, cg in zip(m_rotorset, rotor_cgs)])
        self.final_cg = self.total_m_weighted / self.total_m

    def adjust_cgs_for_ideal(self):
        m_rotorset, rotor_cgs = self.get_rotorset_m_and_cgs()
        total_m = np.sum([self.m_payload, self.m_fuselage,
                          self.m_batteries, self.m_fuelcell,
                          self.m_h2tank, self.wing_1_m, self.wing_2_m] + m_rotorset)
        ideal_total_cg = self.cg_ideal_ac * total_m
        total_m_weighted = np.sum([self.m_payload * self.cg_payload,
                                   self.m_fuselage * self.cg_fuselage,
                                   self.m_batteries * self.cg_batteries,
                                   self.m_fuelcell * self.cg_fuel_cell,
                                   self.m_h2tank * self.cg_h2tank,
                                   self.wing_1_m * self.cg_wing_1,
                                   self.wing_2_m * self.cg_wing_2] + [m * cg for m, cg in zip(m_rotorset, rotor_cgs)])
        adjustment_factor = (ideal_total_cg - total_m_weighted) / total_m

        def adjust_within_bounds(cg, adjustment, min_val, max_val):
            adjusted_cg = cg + adjustment
            return max(min(adjusted_cg, max_val), min_val)

        self.cg_payload = adjust_within_bounds(self.cg_payload, adjustment_factor,  self.payload_cg_min, self.payload_cg_max)
        self.cg_fuselage = adjust_within_bounds(self.cg_fuselage, adjustment_factor, self.fuselage_cg_min, self.fuselage_cg_max)
        self.cg_batteries = adjust_within_bounds(self.cg_batteries, adjustment_factor, self.batteries_cg_min, self.batteries_cg_max)
        self.cg_fuel_cell = adjust_within_bounds(self.cg_fuel_cell, adjustment_factor, self.fuelcell_cg_min, self.fuelcell_cg_max)
        self.cg_h2tank = adjust_within_bounds(self.cg_h2tank, adjustment_factor, self.h2tank_cg_min, self.h2tank_cg_max)
        self.cg_wing_1 = adjust_within_bounds(self.cg_wing_1, adjustment_factor, self.wing_1_cg_min, self.wing_1_cg_max)
        self.cg_wing_2 = adjust_within_bounds(self.cg_wing_2, adjustment_factor, self.wing_2_cg_min, self.wing_2_cg_max)
        self.rotor_cg_positions = [adjust_within_bounds(cg, adjustment_factor, self.rotor_1_cg_min, self.rotor_1_cg_max) if i == 0 else adjust_within_bounds(cg, adjustment_factor, self.rotor_2_cg_min, self.rotor_2_cg_max) for i, cg in enumerate(self.rotor_cg_positions)]
        self.balancing()
        
        # Ensure that h2tank, batteries, and fuelcell CGs are within 30% of each other to avoid having batteries,h2tank, and fuelcell too far apart
        def enforce_component_cg_proximity(self):
            cgs = [self.cg_h2tank, self.cg_batteries, self.cg_fuel_cell]
            min_cg = min(cgs)
            max_cg = max(cgs)
            while (max_cg - min_cg) > 30:
                # Adjust CGs towards their average, within bounds
                avg_cg = sum(cgs) / 3
                self.cg_h2tank = max(min(avg_cg, self.h2tank_cg_max), self.h2tank_cg_min)
                self.cg_batteries = max(min(avg_cg, self.batteries_cg_max), self.batteries_cg_min)
                self.cg_fuel_cell = max(min(avg_cg, self.fuelcell_cg_max), self.fuelcell_cg_min)
                cgs = [self.cg_h2tank, self.cg_batteries, self.cg_fuel_cell]
                min_cg = min(cgs)
                max_cg = max(cgs)
                if (max_cg - min_cg) <= 30:
                    break
                else:
                    break
                # Recalculate final CG
            self.balancing()
        
        while abs(self.final_cg - self.cg_ideal_ac) > 0.01:
            self.adjust_cgs_for_ideal()

    def visual_cg_bar(self):
        components = ['Payload', 'Fuselage', 'Batteries', 'Fuel Cell', 'H2 Tank']
        cg_values = [self.cg_payload, self.cg_fuselage, self.cg_batteries,
                     self.cg_fuel_cell, self.cg_h2tank]
        if self.aircraft_type == 1:
            components += ['Rotors CG']
            cg_values += self.rotor_cg_positions
        elif self.aircraft_type == 2:
            components += [f'Rotor Set {i+1}' for i in range(self.num_rotor_sets)]
            cg_values += self.rotor_cg_positions
        elif self.aircraft_type == 3:
            components += ['Wing & Rotorset 1', 'Wing & Rotorset 2']
            cg_values += [self.cg_wing_1, self.cg_wing_2]
        plt.figure(figsize=(10, 6))
        plt.bar(components, cg_values, color='skyblue', label='System CG Values')
        plt.axhline(y=self.cg_ideal_ac, color='green', linestyle='-', label='Ideal CG = ' + str(self.cg_ideal_ac) + '%')
        plt.axhline(y=self.final_cg, color='red', linestyle='--', label=f'Final CG = {self.final_cg:.2f}%')
        plt.title('System CG Values and Final CG', fontsize=16)
        plt.xlabel('Components', fontsize=14)
        plt.ylabel('CG Position (%)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    def visual_cg_line(self):
        components = ['Payload', 'Fuselage', 'Batteries', 'Fuel Cell', 'H2 Tank']
        cg_values = [self.cg_payload, self.cg_fuselage, self.cg_batteries,
                     self.cg_fuel_cell, self.cg_h2tank]
        masses = [self.m_payload, self.m_fuselage, self.m_batteries,
                  self.m_fuelcell, self.m_h2tank]
        colors = ['blue', 'green', 'orange', 'purple', 'brown']
        
        if self.aircraft_type == 1:
            components += ['Rotors CG']
            cg_values += self.rotor_cg_positions
            masses += [self.rotor_m] * self.num_rotor_sets
            colors += ['grey'] * self.num_rotor_sets
        elif self.aircraft_type == 2:
            components += [f'Rotor Set {i+1}' for i in range(self.num_rotor_sets)]
            cg_values += self.rotor_cg_positions
            masses += [self.rotor_m] * self.num_rotor_sets
            colors += ['grey'] * self.num_rotor_sets
        elif self.aircraft_type == 3:
            components += ['Wing & Rotorset 1', 'Wing & Rotorset 2']
            cg_values += [self.cg_wing_1, self.cg_wing_2]
            masses += [self.wing_1_m, self.wing_2_m]
            colors += ['cyan', 'magenta']
        
        plt.figure(figsize=(12, 6))
        plt.scatter(self.final_cg, 0, color='red', edgecolor='black', s=300, label=f'Final CG = {self.final_cg:.2f}%')
        plt.plot([0, 100], [0, 0], color='black', linewidth=2)
        for i, (cg, mass, color) in enumerate(zip(cg_values, masses, colors)):
            plt.scatter(cg, 0, s=100, color=color, label=components[i])
        plt.yticks([])
        plt.xticks(np.arange(0, 101, 10), fontsize=12)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()

    def input_mass(self):
        self.m_payload = 185  # float(input("Enter payload mass: "))
        self.m_fuselage = 170  # float(input("Enter fuselage mass: "))
        self.m_batteries = 84  # float(input("Enter batteries mass: "))
        self.m_fuelcell = 165  # float(input("Enter fuel cell mass: "))
        self.m_h2tank = 116  # float(input("Enter H2 tank mass: "))
        self.wing_1_m = self.wing_1_m  # float(input("Enter wing 1 mass: "))
        self.wing_2_m = self.wing_1_m  # float(input("Enter wing 2 mass: "))
        self.balancing()
        self.adjust_cgs_for_ideal()

    def print_adjusted_cgs(self):
        print(f"Adjusted CGs:")
        print(f"Payload: {self.cg_payload:.2f}%")
        print(f"Fuselage: {self.cg_fuselage:.2f}%")
        print(f"Batteries: {self.cg_batteries:.2f}%")
        print(f"Fuel Cell: {self.cg_fuel_cell:.2f}%")
        print(f"H2 Tank: {self.cg_h2tank:.2f}%")
        if self.aircraft_type == 3:
            print(f"Wing 1: {self.cg_wing_1:.2f}%")
            print(f"Wing 2: {self.cg_wing_2:.2f}%")
        if self.aircraft_type == 1:
            for i, cg in enumerate(self.rotor_cg_positions):
                print(f"Rotors CG: {self.rotor_cg_positions[0]:.2f}%")
            # print(f"Rotor Set 2: {self.rotor_cg_positions[1]:.2f}%")
        else:
            for i, cg in enumerate(self.rotor_cg_positions):
                print(f"Rotor Set {i+1}: {cg:.2f}%")

if __name__ == '__main__':
    aircraft_type = int(input("Enter aircraft type (1 for Quadrotor, 2 for Coaxial Helicopter, 3 for Compound Aircraft): "))
    MB = MassBalance(aircraft_type)
    MB.input_mass()

    print(f"Final CG: {MB.final_cg:.2f}")

    MB.print_adjusted_cgs()

    MB.visual_cg_bar()
    MB.visual_cg_line()
