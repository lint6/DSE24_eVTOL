import numpy as np
from rotor_sizing import RotorAnalysis  # Import RotorAnalysis class from rotor_sizing

# i am not sure if this is correct
class PerformanceAnalysis:
    def __init__(self, rotor_analysis: RotorAnalysis, vertical_climb_speed: float = 0.76):
        self.rotor_analysis = rotor_analysis
        self.V_vertical_climb = vertical_climb_speed

        # Use values from RotorAnalysis
        self.mtow = rotor_analysis.MTOW_KG
        self.N_rotors = rotor_analysis.N_rotors
        self.rho = rotor_analysis.rho
        self.g = rotor_analysis.g

        # Perform rotor sizing
        self.results = rotor_analysis.perform_analysis(c_t_o_fl=0.1, c_t_o_turn=0.2, c_t_o_turb=0.3)
        self.rotor_radius, self.rotor_diameter, self.v_tip, self.mach_number, self.omega, self.adv_ratio_fl, self.o_fl, self.o_turn, self.o_turb, self.solidity, self.aspect_ratios, self.n_z  = self.results

        # Derived parameters
        self.rotor_area = np.pi * (self.rotor_radius**2)  # Rotor swept area (m^2)
        self.weight_per_rotor = (self.mtow / self.N_rotors) * self.g  # Weight supported by each rotor (N)
        self.FM = 0.78 # figure of merit (assumed)

    def calculate_hover_induced_velocity(self):
        """Calculate hover induced velocity (m/s)."""
        return np.sqrt(self.weight_per_rotor / (2 * self.rho * self.rotor_area))

    def calculate_climb_induced_velocity(self):
        """Calculate climb induced velocity (m/s)."""
        return (-self.V_vertical_climb / 2) + np.sqrt((self.V_vertical_climb / 2)**2 + (self.weight_per_rotor / (2 * self.rho * self.rotor_area)))

    def calculate_disk_area(self):
        """Calculate total disk area of the quadcopter"""
        S_disk = self.rotor_area*4
        return S_disk
    
    def calculate_power_req_hover(self):
        P_req_hover = ((self.mtow*self.g)**(3/2))/(np.sqrt(2*self.rho*self.calculate_disk_area())*self.FM)
        return P_req_hover
    
    def calculate_power_loading_hover(self):
        P_load_hover = (self.mtow) / (self.calculate_power_req_hover()/1000)
        return P_load_hover
    
    
    def display_results(self):
        """Display power requirements and induced velocities."""
        v_i_hover = self.calculate_hover_induced_velocity()
        v_i_climb = self.calculate_climb_induced_velocity()
        disk_area = self.calculate_disk_area()
        power_requirements = self.calculate_power_req_hover()
        power_loading_hover = self.calculate_power_loading_hover()

        #print(f'Hover induced velocity: {v_i_hover:.2f} m/s')
        #print(f'Climb induced velocity: {v_i_climb:.2f} m/s')
        print(f'Disk Area: {disk_area:.2f} m2')
        print(f'Power required for hover: {power_requirements / 1000:.2f} kW')
        print(f'Power loading for hover: {power_loading_hover:.2f} kg/kW')


# Main Execution
if __name__ == "__main__":
    # Initialize RotorAnalysis
    rotor_analysis = RotorAnalysis()

    # Initialize PerformanceAnalysis with the rotor analysis instance
    performance_analysis = PerformanceAnalysis(rotor_analysis)

    # Display results
    performance_analysis.display_results()
