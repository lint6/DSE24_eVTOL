import unittest
from rotor_sizing import RotorAnalysis
import numpy as np

class TestRotorAnalysis(unittest.TestCase):

    #test if initialisation is correct
    def test_default_initialization(self):
        analysis = RotorAnalysis()
        self.assertEqual(analysis.MTOW_KG, 718.89)
        self.assertEqual(analysis.N_rotors, 4)
        self.assertEqual(analysis.rho, 1.225)
        self.assertEqual(analysis.D_v, 0.05*718.89)

    #test if rotor radius calculation is correct for arbitrary mtow
    def test_calculate_radius(self):
        analysis = RotorAnalysis(mtow_kg=1000)
        radius_m = analysis.calculate_radius()
        expected_radius = np.sqrt((1000 * 2.205 / 4) / (3 * np.pi)) / 3.281
        self.assertAlmostEqual(radius_m, expected_radius, places=2)

    #test if tip speed is correctly calculated
    def test_calculate_tip_speed(self):
        analysis = RotorAnalysis()
        tip_speed = analysis.calculate_tip_speed(1.0)  # Test with diameter of 1 meter
        expected_tip_speed = 603 / 3.6  # This is the fixed value

        self.assertEqual(tip_speed, expected_tip_speed)

    #test if mach number is correctly calculated
    def test_calculate_mach_number(self):
        analysis = RotorAnalysis()
        v_max_kmh = 150
        v_tip_ms = 50
        mach_number = analysis.calculate_mach_number(v_max_kmh / 3.6, v_tip_ms)
        expected_mach_number = (v_max_kmh / 3.6 + v_tip_ms) / 343  # Speed of sound = 343 m/s

        self.assertAlmostEqual(mach_number, expected_mach_number, places=3)

    #Test that the perform_analysis function performs calculations without errors and produces expected types of output
    def test_perform_analysis(self):
        analysis = RotorAnalysis()
        # Use realistic thrust coefficients for testing
        results = analysis.perform_analysis(0.1, 0.12, 0.15)

        # Unpack results
        rotor_radius_m, rotor_diameter_m, v_tip_ms, mach_number, omega, adv_ratio_fl, o_fl, o_turn, o_turb, solidity, aspect_ratios, n_z = results

        # Test if results are of expected types and ranges
        self.assertIsInstance(rotor_radius_m, float)
        self.assertIsInstance(rotor_diameter_m, float)
        self.assertIsInstance(v_tip_ms, float)
        self.assertIsInstance(mach_number, float)
        self.assertIsInstance(omega, float)
        self.assertIsInstance(adv_ratio_fl, float)
        self.assertGreater(solidity, 0)
        self.assertTrue(len(aspect_ratios) > 0)

        # Check for positive values
        self.assertGreater(rotor_radius_m, 0)
        self.assertGreater(v_tip_ms, 0)
        self.assertGreater(mach_number, 0)
        self.assertGreater(omega, 0)


    def test_calculate_chord(self):
        analysis = RotorAnalysis()
        # Example with 4 blades, solidity 0.1, and rotor radius of 10 meters
        chord = analysis.calculate_chord(4, 0.1, 10)
        expected_chord = (0.1 * np.pi * 10) / 4

        self.assertAlmostEqual(chord, expected_chord, places=3)

    def test_zero_bank_angle(self):
        analysis = RotorAnalysis(bank_angle=0)
        results = analysis.perform_analysis(0.1, 0.12, 0.15)

        # Check if the load factor in the turn should still be 1 (n_z = 1)
        self.assertEqual(results[11], 1)  # The load factor n_z should be 1


    def test_perform_analysis(self):
        # Create an instance of the analysis class
        analysis = RotorAnalysis()

        # Perform analysis with default coefficients
        c_t_o_fl = 0.1
        c_t_o_turn = 0.2
        c_t_o_turb = 0.3
        results = analysis.perform_analysis(c_t_o_fl, c_t_o_turn, c_t_o_turb)

        # Unpack the results
        rotor_radius_m, rotor_diameter_m, v_tip_ms, mach_number, omega, adv_ratio_fl, o_fl, o_turn, o_turb, solidity, aspect_ratios, n_z = results

        # Test Advance Ratio in forward flight
        omega = v_tip_ms / rotor_radius_m  # Omega calculated manually
        expected_adv_ratio_fl = analysis.V_ne / (omega * rotor_radius_m)
        self.assertAlmostEqual(adv_ratio_fl, expected_adv_ratio_fl, places=2)

        # Test Solidity during Forward Flight
        T_fl = analysis.k_dl * analysis.MTOW_KG * analysis.g  # Thrust in forward flight
        c_t_fl = T_fl / (analysis.N_rotors * (analysis.rho * np.pi * rotor_radius_m**2 * (omega * rotor_radius_m)**2))
        expected_o_fl = c_t_fl / c_t_o_fl
        self.assertAlmostEqual(o_fl, expected_o_fl, places=2)

        # Test Solidity during Turn
        n_z = 1 / np.cos(analysis.bank_angle * (np.pi / 180))  # Load factor
        T_turn = n_z * analysis.k_dl * analysis.MTOW_KG * analysis.g  # Thrust during turn
        c_t_turn = T_turn / (analysis.N_rotors * (analysis.rho * np.pi * rotor_radius_m**2 * (omega * rotor_radius_m)**2))
        expected_o_turn = c_t_turn / c_t_o_turn
        self.assertAlmostEqual(o_turn, expected_o_turn, places=2)

        # Test Solidity during Turbulence
        delta_n = (0.25 * analysis.C_lalpha * (analysis.V_gust / (omega * rotor_radius_m))) / c_t_o_turb
        n_z_turb = 2 + delta_n
        T_turb = n_z_turb * analysis.k_dl * analysis.MTOW_KG * analysis.g  # Thrust during turbulence
        c_t_turb = T_turb / (analysis.N_rotors * (analysis.rho * np.pi * rotor_radius_m ** 2 * (omega * rotor_radius_m) ** 2))
        expected_o_turb = c_t_turb / c_t_o_turb
        self.assertAlmostEqual(o_turb, expected_o_turb, places=2)

        # Test Maximum Solidity
        expected_solidity = max(expected_o_fl, expected_o_turn, expected_o_turb)
        self.assertAlmostEqual(solidity, expected_solidity, places=2)

        # Test Aspect Ratio Calculations (for number of blades 2 to 6)
        for N_bl, chord, AR_bl in aspect_ratios:
            expected_chord = (expected_solidity * np.pi * rotor_radius_m) / N_bl
            self.assertAlmostEqual(chord, expected_chord, places=2)

            expected_AR_bl = rotor_radius_m**2 / (rotor_diameter_m * chord)
            self.assertAlmostEqual(AR_bl, expected_AR_bl, places=2)


if __name__ == "__main__":
    unittest.main()


