import numpy as np
import matplotlib.pyplot as plt

class RotorAnalysis:
    def __init__(self, mtow_kg=718.89, n_rotors=4, dl_imperial=3, bank_angle=30, v_max_kmh=150):
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

    def calculate_radius(self):
        mtow_pounds = self.MTOW_KG * self.POUNDS_PER_KG
        radius_ft = np.sqrt((mtow_pounds / self.N_rotors) / (self.DL_IMPERIAL * np.pi))
        radius_m = radius_ft / self.FT_TO_M
        return radius_m

    def calculate_tip_speed(self, diameter_si):
        return 603 / 3.6  # Tip speed in m/s (fixed)

    def calculate_mach_number(self, v_max, v_tip):
        return (v_max + v_tip) / self.SPEED_OF_SOUND

    def perform_analysis(self, c_t_o_fl, c_t_o_turn, c_t_o_turb):
        rotor_radius_m = self.calculate_radius()
        rotor_diameter_m = 2 * rotor_radius_m

        # Tip speed
        v_tip_ms = self.calculate_tip_speed(rotor_diameter_m)
        mach_number = self.calculate_mach_number(self.V_MAX_KMH / 3.6, v_tip_ms)

        # Rotational speed
        omega = v_tip_ms / rotor_radius_m  # Rotational speed in rad/s

        # Advance ratio
        adv_ratio_fl = self.V_ne / (omega * rotor_radius_m)

        # Thrust and solidity in forward flight
        T_fl = self.k_dl * self.MTOW_KG * self.g
        c_t_fl = T_fl / (self.N_rotors * (self.rho * np.pi * rotor_radius_m**2 * (omega * rotor_radius_m)**2))
        o_fl = c_t_fl / c_t_o_fl

        # Solidity during turn
        n_z = 1 / np.cos(self.bank_angle * (np.pi / 180))  # Load factor
        T_turn = n_z * self.k_dl * self.MTOW_KG * self.g
        c_t_turn = T_turn / (self.N_rotors * (self.rho * np.pi * rotor_radius_m**2 * (omega * rotor_radius_m)**2))
        o_turn = c_t_turn / c_t_o_turn

        # Solidity during turbulence
        delta_n = (0.25 * self.C_lalpha * (self.V_gust / (omega * rotor_radius_m))) / c_t_o_turb
        n_z_turb = 2 + delta_n
        T_turb = n_z_turb * self.k_dl * self.MTOW_KG * self.g
        c_t_turb = T_turb / (self.N_rotors * (self.rho * np.pi * rotor_radius_m ** 2 * (omega * rotor_radius_m) ** 2))
        o_turb = c_t_turb / c_t_o_turb

        # Maximum solidity
        solidity = max(o_fl, o_turn, o_turb)

        # Aspect ratio calculation for different numbers of blades
        N_blades = range(2, 7)
        aspect_ratios = []

        for N_bl in N_blades:
            chord = (solidity * np.pi * rotor_radius_m) / N_bl
            AR_blades = rotor_radius_m**2 / (rotor_diameter_m * chord)
            aspect_ratios.append((N_bl, chord, AR_blades))

        # Return all results including advance ratio
        return rotor_radius_m, rotor_diameter_m, v_tip_ms, mach_number, omega, adv_ratio_fl, o_fl, o_turn, o_turb, solidity, aspect_ratios, n_z

    def plot_aspect_ratios(self, aspect_ratios):
        N_blades = [x[0] for x in aspect_ratios]
        aspect_ratios_values = [x[2] for x in aspect_ratios]

        plt.figure(figsize=(8, 5))
        plt.plot(N_blades, aspect_ratios_values, marker='o', label="Aspect Ratio")
        plt.title("Aspect Ratio vs. Number of Blades")
        plt.xlabel("Number of Blades")
        plt.ylabel("Aspect Ratio")
        plt.grid(True)
        plt.legend()
        plt.show()

    def calculate_chord(self, N_bl, solidity, rotor_radius_m):
        """Calculate the blade chord length given the number of blades and solidity."""
        chord = (solidity * np.pi * rotor_radius_m) / N_bl
        return chord


# Main Execution
if __name__ == "__main__":
    # Initialize the analysis with default values
    analysis = RotorAnalysis()

    # Perform analysis with default coefficients to get advance ratio first
    results = analysis.perform_analysis(c_t_o_fl=0.1, c_t_o_turn=0.2, c_t_o_turb=0.3)

    # Unpack results
    rotor_radius, rotor_diameter, v_tip, mach_number, omega, adv_ratio_fl, o_fl, o_turn, o_turb, solidity, aspect_ratios, n_z = results

    # Print advance ratio and other results
    print(f"Advance Ratio in forward flight = {adv_ratio_fl:.2f}")
    print(f"Rotor Radius [m]: {rotor_radius:.2f}")
    print(f"Rotor Diameter [m]: {rotor_diameter:.2f}")
    print(f"Tip Speed [m/s]: {v_tip:.2f}")
    print(f"Mach Number: {mach_number:.2f}")
    print(f"Rotational Speed (Omega) [rad/s]: {omega:.2f}")
    print(f"RPMs : {omega*(60/(2*np.pi)): .2f}")

    # Now ask for user input for thrust coefficients
    print("Please enter the thrust coefficients:")
    c_t_o_fl = float(input("Enter thrust coefficient for forward flight (c_t_o_fl): "))
    c_t_o_turn = float(input("Enter thrust coefficient for turning flight (c_t_o_turn): "))
    c_t_o_turb = float(input("Enter thrust coefficient for turbulence (c_t_o_turb): "))

    # Perform analysis with user-provided coefficients
    results_with_input = analysis.perform_analysis(c_t_o_fl, c_t_o_turn, c_t_o_turb)

    # Unpack new results
    rotor_radius, rotor_diameter, v_tip, mach_number, omega, adv_ratio_fl, o_fl, o_turn, o_turb, solidity, aspect_ratios, n_z = results_with_input

    # Print key results with user inputs
    print(f"Solidity in Forward Flight: {o_fl:.4f}")
    print(f"Solidity during Turn: {o_turn:.4f}")
    print(f"Solidity during Turbulence: {o_turb:.4f}")
    print(f"Solidity used: {solidity:.4f}")

    # Plot results
    analysis.plot_aspect_ratios(aspect_ratios)

    # Ask user for the number of blades they want and calculate the corresponding chord
    selected_blades = int(input("Enter the number of blades you want (between 2 and 6): "))
    if selected_blades not in range(2, 7):
        print("Invalid number of blades selected. Please choose between 2 and 6 blades.")
    else:
        # Get the selected aspect ratio and chord
        selected_aspect_ratio = None
        selected_chord = None
        for N_bl, chord, AR_bl in aspect_ratios:
            if N_bl == selected_blades:
                selected_aspect_ratio = AR_bl
                selected_chord = chord

        # Print the results
        print(f"Selected number of blades: {selected_blades}")
        print(f"Corresponding blade chord length: {selected_chord:.4f} meters")
        print(f"Corresponding aspect ratio: {selected_aspect_ratio:.4f}")
