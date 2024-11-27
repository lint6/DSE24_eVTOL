import numpy as np
from rotor_sizing import RotorAnalysis  # Import RotorAnalysis class from rotor_sizing

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
        self.rotor_radius, self.rotor_diameter, *_ = self.results

        # Derived parameters
        self.rotor_area = np.pi * (self.rotor_radius**2)  # Rotor swept area (m^2)
        self.weight_per_rotor = (self.mtow / self.N_rotors) * self.g  # Weight supported by each rotor (N)

    def calculate_hover_induced_velocity(self):
        """Calculate hover induced velocity (m/s)."""
        return np.sqrt(self.weight_per_rotor / (2 * self.rho * self.rotor_area))

    def calculate_climb_induced_velocity(self):
        """Calculate climb induced velocity (m/s)."""
        return (-self.V_vertical_climb / 2) + np.sqrt((self.V_vertical_climb / 2)**2 + 
(self.weight_per_rotor / (2 * self.rho * self.rotor_area)))

    def calculate_power_requirements(self):
        """Calculate power requirements for hover and climb."""
        # Induced velocities
        v_i_hover = self.calculate_hover_induced_velocity()
        v_i_climb = self.calculate_climb_induced_velocity()

        # Power requirements per rotor
        P_req_hover_per_rotor = self.weight_per_rotor * v_i_hover
        P_req_climb_per_rotor = self.weight_per_rotor * v_i_climb

        # Total power requirements
        P_req_hover = P_req_hover_per_rotor * self.N_rotors
        P_req_climb = P_req_climb_per_rotor * self.N_rotors

        # Combined hover and climb power
        P_req_hover_climb = (P_req_hover_per_rotor + P_req_climb_per_rotor) * self.N_rotors

        return {
            "hover_power": P_req_hover,
            "climb_power": P_req_climb,
            "combined_power": P_req_hover_climb
        }

    def display_results(self):
        """Display power requirements and induced velocities."""
        v_i_hover = self.calculate_hover_induced_velocity()
        v_i_climb = self.calculate_climb_induced_velocity()
        power_requirements = self.calculate_power_requirements()

        print(f'Hover induced velocity: {v_i_hover:.2f} m/s')
        print(f'Climb induced velocity: {v_i_climb:.2f} m/s')
        print(f'Power required for hover: {power_requirements["hover_power"] / 1000:.2f} kW')
        print(f'Power required for climb: {power_requirements["climb_power"] / 1000:.2f} kW')
        print(f'Power required for hover and climb: {power_requirements["combined_power"] / 1000:.2f} kW')


# Main Execution
if __name__ == "__main__":
    # Initialize RotorAnalysis
    rotor_analysis = RotorAnalysis()

    # Initialize PerformanceAnalysis with the rotor analysis instance
    performance_analysis = PerformanceAnalysis(rotor_analysis)

    # Display results
    performance_analysis.display_results()
