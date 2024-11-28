import numpy as np 
from math import *
import scipy 
import matplotlib.pyplot as plt 


#def calc_p_idd_kappa(self):
#    self.calc_power_ideal_hover()
#    self.calc_kappa_d()
#    self.p_idd_kappa = self.p_idh * self.kappa_d
#    return self.p_idd_kappa

def find_roots(coefficients):
    """
    Finds the roots of a 4th order polynomial equation.

    Parameters:
        coefficients (list): List of coefficients [a, b, c, d, e]
                             for the polynomial ax^4 + bx^3 + cx^2 + dx + e = 0.

    Returns:
        list: Roots of the polynomial.
    """
    # Ensure the coefficients are valid
    if len(coefficients) != 5:
        raise ValueError("You must provide exactly 5 coefficients for a 4th order polynomial.")
    
    # Use numpy's roots function to find the roots
    roots = np.roots(coefficients)
    return roots



class Ducted_Fan_2:
    def __init__(self, mtow, radius=0.625, T_W_R= 1 , V_c=None, density=1.225, D_h0 = 5, k_v_f = 1, V = 100, related_fan = None): 
        # Required Inputs
        self.mtow = mtow  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.D_h0 = D_h0
        self.k_v_f = k_v_f 
        self.V = V

        # Optional Inputs with Defaults
        self.T_W_R = T_W_R  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        
        self.T = None 
        self.mto_weight = None 
        self.rotor_alpha = None
        self.roots = None 

        self.related_fan = related_fan




    def calc_weight(self):
        self.mto_weight = self.mtow * 9.81
        return self.mto_weight
    
    def calc_T(self):
        self.calc_weight()
        self.T = self.mto_weight * np.sqrt(self.k_v_f**2 + (self.D_h0 / self.mto_weight)**2 )
        return self.T
    

    def calc_rotor_alpha(self):
        self.calc_T()
        self.rotor_alpha = - np.arctan(int(self.D_h0) / int(self.T))
        return self.rotor_alpha
    
    def calc_v_f(self):
        self.calc_rotor_alpha()
        my_coefficient = [1, float(-2 * self.V * np.sin(self.rotor_alpha)), (self.V **2), 0, -1 ]
        self.roots = find_roots(coefficients=my_coefficient) * self.related_fan.v_h
        #write an algorithm to chosee THE root you want
        return self.roots 



class Ducted_Fan_1:
    def __init__(self, mtow, radius=0.625, T_W_R= 1 , V_c=None, density=1.225, P_a =410000): 
        # Required Inputs
        self.mtow = mtow  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)

        # Optional Inputs with Defaults
        self.T_W_R = T_W_R  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        self.P_a = P_a

        # Calculated Attributes (initialized as None)
        self.disc_loading = None
        self.thrust_loading = None 
        self.v_h = None  # Hover induced velocity
        self.thrust = None
        #self.hover_induced_velocity = None
        self.v_d = None
        self.v_d_nd = None
        self.v_c = None
        self.v_c_nd = None
        self.kappa_c = None
        self.p_idc = None
        self.p_idh = None 
        self.kappa_p = None 
        self.V_c_kappa = None 
        self.p_idd = None 
        self.kappa_d = None 
        self.p_idd_kappa = None 
        self.V_d_kappa_d = None 


        

    def calc_mtow(self, new_mtow):
        if new_mtow > self.mtow:
            self.mtow = new_mtow
        return self.mtow
    
    #disc loading here is ACTUALLY thurst to area loading, T_W_R is thrust to wright ratio for example 1 is necessary for hovering but not enough to climb
    def calc_disc_loading(self):
        self.disc_loading = (self.mtow * self.T_W_R) / (np.pi * self.radius**2 )
        return self.disc_loading

    def calc_thrust_loading(self):
        self.calc_disc_loading()
        self.thrust_loading = self.disc_loading * 9.81
        return self.thrust_loading

    # Induved velocity in book == v_h 
    def calc_hover_induced_velovity(self):
        self.calc_thrust_loading()
        self.hover_induced_velocity = np.sqrt(self.thrust_loading/(2 * self.density))
        self.v_h = self.hover_induced_velocity
        return self.v_h

    def calc_thrust_hover(self):
        self.calc_hover_induced_velovity()
        self.thrust = (2 * (np.pi * self.radius**2 )) * self.density * (self.v_h)**2
        return self.thrust
    
    # v_d is qquivalnet of v_c but direction is opposite 
    # V_d, V_c is a freestream velocity decent n climb
    # v_i induced velocity (non-hover)
    # has dimension
    # If we are decending fast, V_D/v_d >> 2 , fast_Decent 
    def calc_v_d (self):
        self.calc_disc_loading()
        self.v_d= - 0.5 *self.V_c_kappa - np.sqrt(0.25 * self.V_c_kappa**2 - self.v_h **2)
        return self.v_d
    
    def non_dimensionalize_v_d(self):
        self.calc_v_d()
        self.v_d_nd = (self.v_d/self.v_h)
        return self.v_d_nd
    
    """
    From here we are doing the other side of the graph 
    """
    def calc_v_c (self):
        self.calc_disc_loading()
        self.calc_V_c_kappa()
        self.v_c = - 0.5 * self.V_c_kappa + np.sqrt(0.25 * self.V_c_kappa**2 + self.v_h **2)
        return self.v_c
    
    def non_dimensionalize_v_c(self):
        self.calc_v_c()
        self.v_c_nd = (self.v_c/self.v_h)
        return self.v_c_nd
    
    def calc_kappa_c(self):
        self.calc_disc_loading()
        self.calc_V_c_kappa()
        kappa_v_c = 0.5 * self.V_c_kappa + np.sqrt(0.25 * self.V_c_kappa**2 + self.v_h **2)
        self.kappa_c = kappa_v_c / self.v_h
        return self.kappa_c

    def calc_power_ideal_steady_climb(self):
       self.calc_kappa_c()
       self.p_idc = (self.mtow * 9.81) * self.kappa_c * self.v_h
       return self.p_idc
    
    def calc_power_ideal_hover(self):
        self.calc_thrust_loading()
        self.p_idh = (self.mtow * 9.81) * 1 * self.v_h
        return self.p_idh   
    
    def calc_kappa_p(self):
        self.calc_power_ideal_hover()
        self.kappa_p = self.P_a / self.p_idh
        return self.kappa_p
    
    def calc_V_c_kappa(self):
        self.calc_kappa_p()
        V_c_kappa_nd = self.kappa_p - (1/self.kappa_p)
        self.V_c_kappa = V_c_kappa_nd * self.v_h
        return self.V_c_kappa


    def calc_kappa_d(self):
        self.calc_V_c_kappa()
        self.calc_v_d()
        if np.abs(self.V_c_kappa/self.v_d)>=2 :
            self.kappa_d = 0.5 * (self.V_c_kappa) + np.sqrt(0.25 * self.V_c_kappa**2 - 1)
        elif np.abs(self.V_c_kappa/self.v_d)<1 :
            self.kappa_d = 0.5 * (self.V_c_kappa) + np.sqrt(0.25 * self.V_c_kappa**2 + 1)
        return self.kappa_d


    def calc_p_idd(self):
        self.calc_V_c_kappa()
        self.calc_v_d()
        if np.abs(self.V_c_kappa/self.v_d)>=2 :
            self.p_idd = - self.mtow * 9.81 * (self.v_d + self.V_c_kappa)
        elif np.abs(self.V_c_kappa/self.v_d)<1 :
            self.p_idd = self.mtow * 9.81 * (self.v_d + self.V_c_kappa)
        return self.p_idd

    def calc_V_d_kappa_d(self):
        self.calc_V_c_kappa()
        self.calc_v_d()
        if np.abs(self.V_c_kappa/self.v_d)>=2 :
            V_d_kappa_nd = - (self.kappa_d + (1/self.kappa_d))
            self.V_d_kappa_d = V_d_kappa_nd * self.v_h  
        elif np.abs(self.V_c_kappa/self.v_d)<1 :
            V_d_kappa_nd = - self.kappa_d + (1/self.kappa_d)
            self.V_d_kappa_d = V_d_kappa_nd * self.v_h
        return self.V_d_kappa_d

    

    
fan_1 = Ducted_Fan_1(mtow=float(718/4))

# Calculate hover induced velocity
print(f"Hover Induced Velocity: {fan_1.calc_hover_induced_velovity()}")

# Calculate thrust for hover
print(f"Thrust for Hover: {fan_1.calc_thrust_hover()}")

# Calculate power for steady climb
print(f"Power for Steady Climb: {fan_1.calc_power_ideal_steady_climb()}")

print(f"Rate of Climb rate: {fan_1.calc_V_c_kappa()}")

print(f"calc_power_ideal_hover: {fan_1.calc_power_ideal_hover()}")

print(f"calc_disc_loading: {fan_1.calc_disc_loading()}")

print(f"calc_pidd: {fan_1.calc_p_idd()}")

fan_2 = Ducted_Fan_2(mtow=float(718/4), related_fan=fan_1)
print(f"v from fan 2 from fan 1: {fan_2.calc_v_f()}")


