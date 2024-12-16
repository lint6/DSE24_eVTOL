import unittest
import math
import numpy as np
from unittest.mock import MagicMock
from The_Main import RotorSizing, PowerAnalysis, SoundAnalysis, EnergyAnalysis  # Adjust import based on your file structure

class TestRotorSizing(unittest.TestCase):

    def setUp(self):
        """Set up the test case with default parameters"""
        self.rotor = RotorSizing()

    def test_initialization(self):
        """Test the initialization of the class and check default values"""
        self.assertEqual(self.rotor.MTOW, 718.89)
        self.assertEqual(self.rotor.n_blades, 4)
        self.assertEqual(self.rotor.N_rotors, 4)
        self.assertAlmostEqual(self.rotor.g, 9.80665)
        self.assertEqual(self.rotor.rho, 1.225)
        self.assertEqual(self.rotor.temperature, 15)

    def test_rotor_radius_calculation(self):
        """Test the rotor radius calculation"""
        expected_rotor_radius = np.sqrt((self.rotor.MTOW / self.rotor.N_rotors) / (np.pi * self.rotor.disc_loading))
        self.assertAlmostEqual(self.rotor.rotor_radius, expected_rotor_radius, places=2)

    def test_rotor_diameter_calculation(self):
        """Test the rotor diameter calculation"""
        expected_rotor_diameter = 2 * self.rotor.rotor_radius
        self.assertAlmostEqual(self.rotor.rotor_diameter, expected_rotor_diameter, places=2)

    def test_tip_speed_calculation(self):
        """Test the tip speed calculation"""
        expected_tip_speed = 550 * self.rotor.ft_to_m
        self.assertAlmostEqual(self.rotor.tip_speed, expected_tip_speed, places=2)

    def test_omega_calculation(self):
        """Test the angular velocity (omega) calculation"""
        expected_omega = self.rotor.tip_speed / self.rotor.rotor_radius
        self.assertAlmostEqual(self.rotor.omega, expected_omega, places=2)

    def test_RPM_calculation(self):
        """Test the RPM calculation"""
        expected_RPM = self.rotor.omega * (60 / 2 * np.pi)
        self.assertAlmostEqual(self.rotor.RPM, expected_RPM, places=2)

    def test_max_forward_velocity(self):
        """Test the max forward velocity calculation"""
        expected_max_forward_velocity = (self.rotor.max_tip_mach * self.rotor.speed_of_sound) - self.rotor.tip_speed
        self.assertAlmostEqual(self.rotor.max_forward_velocity, expected_max_forward_velocity, places=2)

    def test_maximum_solidity(self):
        """Test the maximum solidity calculation"""
        self.rotor.update_parameters()  # Ensure all parameters are updated
        expected_maximum_solidity = np.max([self.rotor.solidity_forward_flight, self.rotor.solidity_turn, self.rotor.solidity_turbulent])
        self.assertAlmostEqual(self.rotor.maximum_solidity, expected_maximum_solidity, places=2)

    def test_chord_calculation(self):
        """Test the blade chord calculation"""
        self.rotor.update_parameters()  # Ensure all parameters are updated
        expected_chord = (self.rotor.maximum_solidity * np.pi * self.rotor.rotor_radius) / (self.rotor.coaxial * self.rotor.n_blades)
        self.assertAlmostEqual(self.rotor.chord, expected_chord, places=3)

    def test_aspect_ratio_calculation(self):
        """Test the aspect ratio calculation"""
        self.rotor.update_parameters()  # Ensure all parameters are updated
        expected_aspect_ratio = self.rotor.rotor_radius**2 / (self.rotor.rotor_radius * self.rotor.chord)
        self.assertAlmostEqual(self.rotor.aspect_ratio, expected_aspect_ratio, places=2)

    def test_iterate_design(self):
        """Test the iterate_design method"""
        initial_chord = self.rotor.chord
        initial_aspect_ratio = self.rotor.aspect_ratio
        self.rotor.iterate_design(new_MTOW=800, new_n_blades=6)  # Update design with new parameters
        self.assertNotEqual(self.rotor.chord, initial_chord)
        self.assertNotEqual(self.rotor.aspect_ratio, initial_aspect_ratio)

    def test_visual_blade_vs_aspect_ratio(self):
        """Test the visual_blade_vs_aspect_ratio method"""
        # Since this method plots a graph, we can't assert the visual output directly,
        # but we can check if the plot is generated without errors
        try:
            self.rotor.visual_blade_vs_aspect_ratio()
        except Exception as e:
            self.fail(f"visual_blade_vs_aspect_ratio() raised {e} unexpectedly")


class TestPowerAnalysis(unittest.TestCase):

    def setUp(self):
        # Create a mock RotorSizing object to pass to PowerAnalysis
        self.mock_rotorsizing = RotorSizing()
        self.mock_rotorsizing.FM = 0.1
        self.mock_rotorsizing.maximum_solidity = 0.1
        self.mock_rotorsizing.rho = 1.225  # Air density at sea level in kg/m^3
        self.mock_rotorsizing.omega = 300  # Rotor speed in rad/s
        self.mock_rotorsizing.rotor_radius = 2.5  # Rotor radius in meters
        self.mock_rotorsizing.MTOW = 1500  # Maximum Takeoff Weight in kg
        self.mock_rotorsizing.k_dl = 1.2  # Downwash factor
        self.mock_rotorsizing.k_int = 1.28  # Induced power correction factor
        self.mock_rotorsizing.lift_slope = 5.0  # Lift slope for NACA 0012 airfoil
        self.mock_rotorsizing.A_eq = 0.5  # Equivalent flat plate area in m^2
        self.mock_rotorsizing.T_forward_flight = 500  # Thrust for forward flight in N
        self.mock_rotorsizing.N_rotors = 4  # Number of rotors

        # Instantiate PowerAnalysis with the mock object
        self.pa = PowerAnalysis(self.mock_rotorsizing)

    def test_initialization(self):
        # Test if the initialization sets parameters correctly
        self.assertEqual(self.pa.g, 9.80665)  # Gravitational acceleration
        self.assertEqual(self.pa.MTOW, 1500)  # Maximum Takeoff Weight
        self.assertEqual(self.pa.rho, 1.225)  # Air density at sea level
        self.assertEqual(self.pa.omega, 300)  # Rotor speed in rad/s
        self.assertEqual(self.pa.rotor_radius, 2.5)  # Rotor radius
        self.assertEqual(self.pa.V_point, 13)  # Cruise speed

    def test_forward_flight(self):
        # Test the forward flight calculations
        self.pa.forward_flight()

        # Check the profile drag power (P_p)
        self.assertIsNotNone(self.pa.P_p)
        self.assertGreater(self.pa.P_p, 0)

        # Check the induced power (P_i)
        self.assertIsNotNone(self.pa.P_i)
        self.assertGreater(self.pa.P_i, 0)

        # Check the parasitic power (P_par)
        self.assertIsNotNone(self.pa.P_par)
        self.assertGreater(self.pa.P_par, 0)

        # Check the total power for level flight (P_total_level)
        self.assertIsNotNone(self.pa.P_total_level)
        self.assertGreater(self.pa.P_total_level, 0)

    def test_min_power(self):
        # Test the calculation of minimum power
        self.pa.plot_power_components()

        self.assertIsNotNone(self.pa.min_power)
        self.assertGreater(self.pa.min_power, 0)
        self.assertIsNotNone(self.pa.min_power_velocity)

    def test_iterate_design(self):
        # Test that the iterate_design method updates the power calculations correctly
        new_MTOW_N = 18000  # New weight in N
        new_V_point = 15  # New cruise speed

        # Run the iteration
        self.pa.iterate_design(new_MTOW_N=new_MTOW_N, new_V_point=new_V_point)

        # Check that the new parameters are correctly updated
        self.assertEqual(self.pa.MTOW_N, new_MTOW_N)
        self.assertEqual(self.pa.V_point, new_V_point)

    def test_power_values(self):
        # Test if power values for climb and descent are computed
        self.pa.forward_flight()

        # Check vertical climb power
        self.assertIsNotNone(self.pa.P_vertical_climb)
        self.assertGreater(self.pa.P_vertical_climb, 0)

        # Check vertical descent power
        self.assertIsNotNone(self.pa.P_vertical_descent)
        self.assertGreater(self.pa.P_vertical_descent, 0)

        # Check HOGE power
        self.assertIsNotNone(self.pa.P_hoge)
        self.assertGreater(self.pa.P_hoge, 0)

    def test_plot(self):
        # Test the plotting function (check if it doesn't throw errors)
        try:
            self.pa.plot_power_components()
        except Exception as e:
            self.fail(f"Plotting failed with exception: {e}")


class TestSoundAnalysis(unittest.TestCase):

    def setUp(self):
        # Create a mock RotorSizing object to pass to SoundAnalysis
        self.mock_rotorsizing = RotorSizing()
        self.mock_rotorsizing.rotor_radius = 10  # Rotor radius in meters
        self.mock_rotorsizing.omega = 300  # Rotor rotational speed in rad/s
        self.mock_rotorsizing.V_max = 150  # Max flight speed in m/s
        self.mock_rotorsizing.speed_of_sound = 343  # Speed of sound in m/s
        self.mock_rotorsizing.n_blades = 3  # Number of blades per rotor
        self.mock_rotorsizing.T_forward_flight = 1000  # Thrust in Newtons
        self.mock_rotorsizing.N_rotors = 4  # Number of rotors
        self.mock_rotorsizing.coaxial = 1  # Coaxial rotors flag

        # Instantiate SoundAnalysis with the mock object
        self.sa = SoundAnalysis(self.mock_rotorsizing)

    def test_initialization(self):
        # Test if the initialization sets parameters correctly
        self.assertEqual(self.sa.m_to_f, 1/0.3048)  # Conversion factor for meters to feet
        self.assertEqual(self.sa.N_to_lbs, (1/9.80665)*2.20462)  # Conversion factor for Newtons to pounds
        self.assertEqual(self.sa.R, self.mock_rotorsizing.rotor_radius * self.sa.m_to_f)  # Rotor radius in feet
        self.assertEqual(self.sa.A, np.pi * (self.sa.R ** 2))  # Rotor area in ft^2
        self.assertEqual(self.sa.T, self.mock_rotorsizing.T_forward_flight * self.sa.N_to_lbs / self.mock_rotorsizing.N_rotors)  # Thrust in lbs

    def test_rotational_noise(self):
        # Test the calculation of rotational noise SPL
        self.sa.rotational_noise()

        # Check the calculated values
        self.assertIsNotNone(self.sa.rotational_SPL)
        self.assertGreater(self.sa.rotational_SPL, 0)
        self.assertIsNotNone(self.sa.rotational_SPL_total)
        self.assertGreater(self.sa.rotational_SPL_total, 0)
        self.assertIsNotNone(self.sa.f_rotational)
        self.assertGreater(self.sa.f_rotational, 0)

    def test_vortex_noise(self):
        # Test the calculation of vortex noise SPL
        self.sa.vortex_noise()

        # Check the calculated values
        self.assertIsNotNone(self.sa.vortex_SPL)
        self.assertGreater(self.sa.vortex_SPL, 0)
        self.assertIsNotNone(self.sa.vortex_SPL_total)
        self.assertGreater(self.sa.vortex_SPL_total, 0)
        self.assertIsNotNone(self.sa.f_vortex)
        self.assertGreater(self.sa.f_vortex, 0)

    def test_display_parameters_rotor(self):
        # Test display function for rotor-related parameters
        # We'll capture printed output and check if the output contains expected parameters
        from io import StringIO
        import sys

        captured_output = StringIO()
        sys.stdout = captured_output
        self.sa.display_parameters_rotor()
        sys.stdout = sys.__stdout__

        # Check if output contains expected parameters
        self.assertIn("This configuration has 4 rotors", captured_output.getvalue())
        self.assertIn("Single rotor SPL:", captured_output.getvalue())
        self.assertIn("All rotor SPL:", captured_output.getvalue())
        self.assertIn("Fundamental frequency:", captured_output.getvalue())

    def test_display_parameters_vortex(self):
        # Test display function for vortex-related parameters
        # We'll capture printed output and check if the output contains expected parameters
        from io import StringIO
        import sys

        captured_output = StringIO()
        sys.stdout = captured_output
        self.sa.display_paramenters_vortex()
        sys.stdout = sys.__stdout__

        # Check if output contains expected parameters
        self.assertIn("Single rotor SPL:", captured_output.getvalue())
        self.assertIn("All rotor SPL:", captured_output.getvalue())
        self.assertIn("Vortex frequency:", captured_output.getvalue())


class TestEnergyAnalysis(unittest.TestCase):

    def setUp(self):
        # Create a mock PowerAnalysis instance
        self.mock_power_analysis = MagicMock()
        self.mock_power_analysis.P = {
            'HIGE1': 100, 'V_climb': 150, 'HOGE1': 200, 'Climb1': 300,
            'Cruise1': 400, 'Descent1': 250, 'HOGE2': 200, 'Loiter': 100,
            'Climb2': 300, 'Cruise2': 350, 'Descent2': 250, 'HOGE3': 200,
            'V_Descent': 150, 'HIGE2': 100
        }
        self.mock_power_analysis.rho = 1.225
        self.mock_power_analysis.vertical_climb = 5.0
        self.mock_power_analysis.min_power_velocity = 50
        self.mock_power_analysis.min_power_velocity_CD = 60
        self.mock_power_analysis.min_power_velocity_descent = 40
        self.mock_power_analysis.gamma_CD = 10
        self.mock_power_analysis.gamma_descent = 5
        
        # Create an instance of EnergyAnalysis with the mock PowerAnalysis
        self.energy_analysis = EnergyAnalysis(power=self.mock_power_analysis)

    def test_calculate_missionphase_time(self):
        # Run the time calculation
        times = self.energy_analysis.calculate_missionphase_time()
        
        # Test if the times for specific phases are correctly calculated
        self.assertEqual(times['V_climb'], 60.0 / self.mock_power_analysis.vertical_climb)  # Vertical climb time
        self.assertEqual(times['Climb1'], (300.0 - 60.0) / np.tan(np.radians(10)) / self.mock_power_analysis.min_power_velocity_CD)
        self.assertEqual(times['Descent1'], (300.0 - 60.0) / 7.6)  # Descent time
        self.assertEqual(times['Cruise1'], 30000 / self.mock_power_analysis.min_power_velocity)  # Cruise 1 time

    def test_calculate_energy_required(self):
        # Run the energy calculation
        energies = self.energy_analysis.calculate_energy_required()

        # Check that energy for a phase is correctly calculated
        self.assertAlmostEqual(energies['V_climb'], (self.energy_analysis.mission_data['times']['V_climb'] / 3600) * 150 / 0.85)
        self.assertAlmostEqual(energies['Climb1'], (self.energy_analysis.mission_data['times']['Climb1'] / 3600) * 300 / 0.85)

    def test_calculate_amps(self):
        # Run the amperage calculation
        amps = self.energy_analysis.calculate_amps()

        # Test that amps are calculated for a phase
        self.assertAlmostEqual(amps['V_climb'], 150 / (840 * 0.85))  # V_climb amperage

        # Check if max amperage is correctly calculated
        self.assertEqual(amps['max'], max(amps.values()))

    def test_visual_PEMFC_power(self):
        # This function doesn't return anything, but we can check if it runs without error
        try:
            self.energy_analysis.visual_PEMFC_power()
            success = True
        except Exception as e:
            success = False
            print(f"Error: {e}")
        
        self.assertTrue(success)

    def test_visual_PEMFC_energy(self):
        # Similar to the previous function, we just check if it runs without error
        try:
            self.energy_analysis.visual_PEMFC_energy()
            success = True
        except Exception as e:
            success = False
            print(f"Error: {e}")
        
        self.assertTrue(success)


if __name__ == "__main__":
    unittest.main()
