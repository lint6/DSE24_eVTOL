# Made by Ibrahim and Floris

import numpy as np
import matplotlib.pyplot as plt

class MassBalance:

    def __init__(self, aircraft_type):
        self.aircraft_type = aircraft_type
        self.MTOW = 718.89
        self.m_payload = 185
        self.m_airframe = 152
        self.m_batteries = 84
        self.m_fuelcell = 165
        self.m_h2tank = 114
        if self.aircraft_type == 1:
            #Define Ideal Total CG for Quadrotor
            self.cg_ideal_ac = 50
            # Define ideal min/max CG positions for each component
            self.payload_cg_min = 10
            self.payload_cg_max = 20
            self.airframe_cg_min = 45
            self.airframe_cg_max = 55
            self.batteries_cg_min = 30
            self.batteries_cg_max = 70
            self.fuelcell_cg_min = 50
            self.fuelcell_cg_max = 90
            self.h2tank_cg_min = 35
            self.h2tank_cg_max = 70
            # Define rotor set parameters
            self.num_rotor_sets = 2
            self.rotor_m = 5
            self.rotor_1_cg_min = 20
            self.rotor_1_cg_max = 30
            self.rotor_2_cg_min = 70
            self.rotor_2_cg_max = 80
            self.rotor_cg_positions = [self.rotor_1_cg_min, self.rotor_2_cg_max]
        elif self.aircraft_type == 2:
            # Define parameters for Coaxial Helicopter
            self.cg_ideal_ac = 39 
            self.payload_cg_min = 11
            self.payload_cg_max = 21
            self.airframe_cg_min = 30
            self.airframe_cg_max = 45
            self.batteries_cg_min = 31
            self.batteries_cg_max = 56
            self.fuelcell_cg_min = 31
            self.fuelcell_cg_max = 56
            self.h2tank_cg_min = 31
            self.h2tank_cg_max = 56
            self.num_rotor_sets = 1
            self.rotor_m = 36 #Based on 3.9 m radius from this site https://aircommand.com/pages/rotor-blade-selection-and-planning#
            self.rotor_1_cg_min = self.cg_ideal_ac 
            self.rotor_1_cg_max = self.cg_ideal_ac 
            self.rotor_cg_positions = [self.cg_ideal_ac]
            pass
        elif self.aircraft_type == 3:
            # Define parameters for Compound Aircraft
            self.cg_ideal_ac = 50 #ASK COMPOUND GROUP
            self.num_rotor_sets = 2
            self.rotor_m = 5 #ASK COMPOUND GROUP
            self.rotor_1_cg_min = 20 #ASK COMPOUND GROUP
            self.rotor_1_cg_max = 30 #ASK COMPOUND GROUP
            self.rotor_2_cg_min = 70
            self.rotor_2_cg_max = 80
            self.rotor_cg_positions = [self.rotor_1_cg_min, self.rotor_2_cg_max]
            pass
        else:
            raise ValueError("Invalid aircraft type. Choose 1 for Quadrotor, 2 for Coaxial Helicopter, or 3 for Compound Aircraft.")

        self.cg_payload = np.random.uniform(10, 20)
        self.cg_airframe = np.random.uniform(45, 55)
        self.cg_batteries = np.random.uniform(30, 70)
        self.cg_fuel_cell = np.random.uniform(50, 90)
        self.cg_h2tank = np.random.uniform(35, 70)

        self.balancing()

    def get_rotorset_m_and_cgs(self):
        m_rotorset = [self.rotor_m] * self.num_rotor_sets * 2
        cgs = []
        for cg in self.rotor_cg_positions:
            cgs.extend([cg] * self.num_rotor_sets)
        return m_rotorset, cgs

    def balancing(self):
        m_rotorset, rotor_cgs = self.get_rotorset_m_and_cgs()
        self.total_m = np.sum([self.m_payload, self.m_airframe,
                               self.m_batteries, self.m_fuelcell,
                               self.m_h2tank] + m_rotorset)
        self.total_m_weighted = np.sum([self.m_payload * self.cg_payload,
                                        self.m_airframe * self.cg_airframe,
                                        self.m_batteries * self.cg_batteries,
                                        self.m_fuelcell * self.cg_fuel_cell,
                                        self.m_h2tank * self.cg_h2tank] + [m * cg for m, cg in zip(m_rotorset, rotor_cgs)])
        self.final_cg = self.total_m_weighted / self.total_m

    def adjust_cgs_for_ideal(self):
        m_rotorset, rotor_cgs = self.get_rotorset_m_and_cgs()
        total_m = np.sum([self.m_payload, self.m_airframe,
                          self.m_batteries, self.m_fuelcell,
                          self.m_h2tank] + m_rotorset)
        ideal_total_cg = self.cg_ideal_ac * total_m
        total_m_weighted = np.sum([self.m_payload * self.cg_payload,
                                   self.m_airframe * self.cg_airframe,
                                   self.m_batteries * self.cg_batteries,
                                   self.m_fuelcell * self.cg_fuel_cell,
                                   self.m_h2tank * self.cg_h2tank] + [m * cg for m, cg in zip(m_rotorset, rotor_cgs)])
        adjustment_factor = (ideal_total_cg - total_m_weighted) / total_m

        def adjust_within_bounds(cg, adjustment, min_val, max_val):
            adjusted_cg = cg + adjustment
            return max(min(adjusted_cg, max_val), min_val)

        self.cg_payload = adjust_within_bounds(self.cg_payload, adjustment_factor,  self.payload_cg_min, self.payload_cg_max)
        self.cg_airframe = adjust_within_bounds(self.cg_airframe, adjustment_factor, self.airframe_cg_min, self.airframe_cg_max)
        self.cg_batteries = adjust_within_bounds(self.cg_batteries, adjustment_factor, self.batteries_cg_min, self.batteries_cg_max)
        self.cg_fuel_cell = adjust_within_bounds(self.cg_fuel_cell, adjustment_factor, self.fuelcell_cg_min, self.fuelcell_cg_max)
        self.cg_h2tank = adjust_within_bounds(self.cg_h2tank, adjustment_factor, self.h2tank_cg_min, self.h2tank_cg_max)
        self.rotor_cg_positions = [adjust_within_bounds(cg, adjustment_factor, self.rotor_1_cg_min, self.rotor_1_cg_max) if i == 0 else adjust_within_bounds(cg, adjustment_factor, self.rotor_2_cg_min, self.rotor_2_cg_max) for i, cg in enumerate(self.rotor_cg_positions)]

        self.balancing()
        
        while abs(self.final_cg - self.cg_ideal_ac) > 0.01:
            self.adjust_cgs_for_ideal()

    def visual_cg_bar(self):
        components = ['Payload', 'Airframe', 'Batteries', 'Fuel Cell', 'H2 Tank'] + [f'Rotor Set {i+1}' for i in range(self.num_rotor_sets)]
        cg_values = [self.cg_payload, self.cg_airframe, self.cg_batteries,
                     self.cg_fuel_cell, self.cg_h2tank] + self.rotor_cg_positions
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
        components = ['Payload', 'Airframe', 'Batteries', 'Fuel Cell', 'H2 Tank'] + [f'Rotor Set {i+1}' for i in range(self.num_rotor_sets)]
        cg_values = [self.cg_payload, self.cg_airframe, self.cg_batteries,
                     self.cg_fuel_cell, self.cg_h2tank] + self.rotor_cg_positions
        masses = [self.m_payload, self.m_airframe, self.m_batteries,
                  self.m_fuelcell, self.m_h2tank] + [self.rotor_m * 2] * self.num_rotor_sets
        plt.figure(figsize=(12, 6))
        plt.plot([0, 100], [0, 0], color='black', linewidth=2)
        plt.scatter(cg_values, np.zeros_like(cg_values), s=np.array(masses) * 2, color='blue', label='Component CG')
        for i, label in enumerate(components):
            plt.text(cg_values[i], 0.05, label, ha='center', fontsize=10)
        plt.scatter(self.final_cg, 0, color='red', s=200, label=f'Final CG = {self.final_cg:.2f}%')
        plt.yticks([])
        plt.xticks(np.arange(0, 101, 10), fontsize=12)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()

    def input_mass(self):
        self.m_payload = 185  # float(input("Enter payload mass: "))
        self.m_airframe = 170  # float(input("Enter airframe mass: "))
        self.m_batteries = 84  # float(input("Enter batteries mass: "))
        self.m_fuelcell = 165  # float(input("Enter fuel cell mass: "))
        self.m_h2tank = 116  # float(input("Enter H2 tank mass: "))
        self.balancing()
        self.adjust_cgs_for_ideal()

    def print_adjusted_cgs(self):
        print(f"Adjusted CGs:")
        print(f"Payload: {self.cg_payload:.2f}%")
        print(f"Airframe: {self.cg_airframe:.2f}%")
        print(f"Batteries: {self.cg_batteries:.2f}%")
        print(f"Fuel Cell: {self.cg_fuel_cell:.2f}%")
        print(f"H2 Tank: {self.cg_h2tank:.2f}%")
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
