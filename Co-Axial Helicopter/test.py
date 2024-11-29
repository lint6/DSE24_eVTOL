import numpy as np
import unittest 
from co_rotor_sizing import RotorSizing

class TestRotorSizing(unittest.TestCase):

    def test_initialization(self):
    # test initialization 
        sizing = RotorSizing()
        self.assertequal(sizing.MTOW, 718.89)
        self.assertequal(sizing.n_blades, 6)
        self.assertequal(sizing.max_tip_mach, 0.85)
    
    def test_update_parameters(self):
        sizing = RotorSizing()
        self.assertIsInstance

        


    

