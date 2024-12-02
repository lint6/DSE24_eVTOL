"""
import numpy as np
import matplotlib.pyplot as plt
from peformance import PerformanceAnalysis  # Import PerformanceAnalysis class from performance
from rotor_sizing import RotorAnalysis  # Import RotorAnalysis class from rotor_sizing

class PowerAnalysis:
    def __init__(self, performance: PerformanceAnalysis, rotor_analysis: RotorAnalysis):
        #values from rotor sizing
        self.N_rotors = rotor_analysis.N_rotors     # get the amount of rotors    
        #values from performance
        self.P_req_hover = performance.calculate_power_req_hover()     # get the required power for hover
        #add power required for take off, climb, cruise, descent, loiter, land
        self.t1 = None 
        self.t2
    
    def calculate_missionphase_time(self):
        t_1 = 15 
        self.t_1 = t_1 
        return self.t_1 
    
    def calcualte_thurst(self, P, t):
        self.calculate_pahere()
        E = P*t

"""
