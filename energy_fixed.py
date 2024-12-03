import numpy as np
from The_Main import RotorSizing, PowerAnalysis

# Placeholder for the Power class (this should be replaced with actual data from performance data).
class Power:
    def __init__(self):
        self.P = {
            'HIGE1': 60000,
            'V_climb': 85000,
            'HOGE1': 65000,
            'Climb1': 100000,
            'Cruise1': 75000,
            'Descent1': 40000,
            'HOGE2': 65000,
            'Loiter': 70000,
            'Climb2': 100000,
            'Cruise2': 75000,
            'Descent2': 40000,
            'HOGE3': 65000,
            'V_Descent': 37500,
            'HIGE2': 60000
        }

class PowerAnalysis:
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


# Assuming the above classes Power and PowerAnalysis are defined in the same script or imported from another module

def main():
    # Create an instance of PowerAnalysis
    power_analysis = PowerAnalysis()

    # Calculate mission phase times
    times = power_analysis.calculate_missionphase_time()
    print("Mission Phase Times (in seconds):")
    for phase, time in times.items():
        print(f"{phase}: {time} seconds")
    print(f"\nTotal Time: {times['total']} seconds\n")

    # Calculate energy requirements
    energies = power_analysis.calculate_energy_required()
    print("Energy Required (in Wh):")
    for phase, energy in energies.items():
        print(f"{phase}: {energy:.2f} Wh")
    print(f"\nTotal Energy: {energies['total']:.2f} Wh\n")

    # Calculate current (amps)
    amps = power_analysis.calculate_amps()
    print("Current (Amps) per Phase:")
    for phase, current in amps.items():
        if phase != 'max':  # Skip 'max' as it's a summary
            print(f"{phase}: {current:.2f} A")
    print(f"\nMaximum Current: {amps['max']:.2f} A\n")


if __name__ == "__main__":
    main()
