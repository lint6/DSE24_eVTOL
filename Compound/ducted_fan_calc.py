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
    if len(coefficients) != 5:
        raise ValueError("You must provide exactly 5 coefficients for a 4th order polynomial.")
    roots = np.roots(coefficients)  
    real_roots = roots[np.isreal(roots)]
    real_roott = min(real_roots[real_roots>0])
    real_root = real_roott.real
    return real_root



class Ducted_Fan_1: #vertical flight
    def __init__(self, mass, radius=0.625, T_W_R= 1 , V_c=None, density=1.225, P_a = 58360): 
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)

        # Optional Inputs with Defaults
        self.T_W_R = T_W_R  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        self.P_a = P_a # Power available set by user

        # Calculated Attributes (initialized as None)
        self.disc_loading = None
        self.thrust_loading = None 
        self.v_h = None  # Hover induced velocity
        self.thrust = None # thurst
        #self.hover_induced_velocity = None
        self.v_d = None #induced velocitt decent
        self.v_d_nd = None #above non dimensional
        self.v_c = None #induced velocity climb
        self.v_c_nd = None #above non dimensional
        self.kappa_c = None # #ratio of ideal power abaliable o ideal power required in hover, from climb rate
        self.p_idc = None # power ideal climb
        self.p_idh = None  # power ideal hover
        self.kappa_p = None  # #ratio of ideal power abaliable o ideal power required in hover, from power avalibale 
        self.V_c_kappa = None # Vertical climb rate from kappa 
        self.p_idd = 34000 # power ideal decent
        self.kappa_d = None  # #ratio of ideal power abaliable o ideal power required in hover, from power avalable
        self.p_idd_kappa = None # power ideal descent from, kappa 
        self.V_d_kappa_d = None  #vertical Descnet rate from kappa RATE


        

    def calc_mass(self, new_mass):
        if new_mass > self.mass:
            self.mass = new_mass
        return self.mass
    
    #disc loading here is ACTUALLY thurst to area loading, T_W_R is thrust to wright ratio for example 1 is necessary for hovering but not enough to climb
    def calc_disc_loading(self):
        self.disc_loading = (self.mass * self.T_W_R) / (np.pi * self.radius**2 )
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
    
    def calc_kappa_c(self): #ratio of ideal power abaliable o ideal power required in hover
        self.calc_disc_loading()
        self.calc_V_c_kappa()
        kappa_v_c = 0.5 * self.V_c_kappa + np.sqrt(0.25 * self.V_c_kappa**2 + self.v_h **2)
        self.kappa_c = kappa_v_c / self.v_h
        return self.kappa_c

    def calc_power_ideal_steady_climb(self):
       self.calc_kappa_c()
       self.p_idc = (self.mass * 9.81) * self.kappa_c * self.v_h
       return self.p_idc
    
    def calc_power_ideal_hover(self):
        self.calc_thrust_loading()
        self.p_idh = (self.mass * 9.81) * 1 * self.v_h
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
            self.p_idd = - self.mass * 9.81 * (self.v_d + self.V_c_kappa)
        elif np.abs(self.V_c_kappa/self.v_d)<1 :
            self.p_idd = self.mass * 9.81 * (self.v_d + self.V_c_kappa)
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

    
class Ducted_Fan_2: #pure horizontal
    def __init__(self, mass, Cd0, V, radius=0.625, T_W_R= 1 , gamma = 0, V_c=None, density=1.225, k_v_f = 1, related_fan = None): 
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.Cd0 = Cd0
        self.k_v_f = k_v_f 
        self.V = V # we need to know this somehow first .... 
        self.gamma = gamma

        # Optional Inputs with Defaults
        self.T_W_R = T_W_R  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        
        self.D = None
        self.D_h0 = None 
        self.V_nd = None
        self.T = None  # Thurst
        self.mto_weight = None #weight
        self.rotor_alpha = None #rotor blade pitch angle
        self.v_f_root = None #induced velocity forward flgiht 
        self.V_hor = None #horizontal velocity.
        self.p_idf = None  #steady horizontal flgith ideal power for rotor

        self.related_fan = related_fan

    def calc_V_horizontal(self):
        self.V_hor = self.V * np.cos(self.gamma)
        return self.V_hor
    
    def calc_V_nd(self):
        self.V_nd = self.V / self.related_fan.v_h
        return self.V_nd
    
    def calc_D(self):
        self.D = (0.5) * self.density * (self.V **2) * self.Cd0 * (self.radius**2 * np.pi)
        self.D_h0 = self.D * np.cos(self.gamma)
        return self.D, self.D_h0

    def calc_weight(self):
        self.mto_weight = self.mass * 9.80665
        return self.mto_weight
    
    def calc_T(self):
        self.calc_D()
        self.calc_weight()
        self.T = self.mto_weight * np.sqrt(self.k_v_f**2 + (self.D_h0 / self.mto_weight)**2 )
        return self.T
    
    def calc_rotor_alpha(self):
        self.calc_T()
        self.calc_D()
        self.rotor_alpha = - np.arctan(int(self.D_h0) / int(self.T))
        # print( "rotor alpha is", self.rotor_alpha*180/np.pi )
        return self.rotor_alpha
    
    def calc_v_f(self):
        self.calc_rotor_alpha()
        self.calc_V_nd()
        my_coefficient = [1, (-2 * self.V_nd * np.sin(self.rotor_alpha)), (self.V_nd **2), 0, -1 ]
        # print("The coefficients are ", my_coefficient)
        self.v_f_root = find_roots(coefficients=my_coefficient) * self.related_fan.v_h
        return self.v_f_root
    
    
    def calc_p_idf(self):
        self.calc_D()
        self.calc_T()
        self.calc_V_horizontal()
        self.calc_rotor_alpha()
        self.calc_v_f()
        
        power_induced = self.T * self.v_f_root
        power_parasitic = self.V_hor * self.D_h0
        self.p_idf = power_induced + power_parasitic
        # print("v_horis ", self.V_hor)
        # print("drag is", self.D_h0)
        # print("thrust is", self.T)
        # print("horizontal induced velocity is", self.v_f_root)
        # print(f'Power: {self.p_idf}')
        # print(f'power_induced: {power_induced}')
        # print(f'power_parasitic: {power_parasitic}')
        # print('----------------------------------')
        return self.p_idf, power_induced, power_parasitic
  


class Ducted_Fan_3: #angled climb
    def __init__(self, mass, gamma, radius=0.625, T_W_R= 1 , V_c=None, V =4, density=1.225, D_h0 = 500, k_v_f = 1, P_a =45000,  related_fan1 = None, related_fan2 = None): 
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.D_h0 = D_h0 #Drag horizontal
        self.k_v_f = k_v_f  #coefficient of vertical drag ~ 1
        self.V = V # Velocity we give
        self.gamma = gamma # airpath angle
        self.P_a = P_a # power avilable

        # Optional Inputs with Defaults
        self.T_W_R = T_W_R  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        
        #self.T = fan_2.T
        #self.mto_weight = fan_2.mto_weight
        #self.rotor_alpha = fan_2.rotor_alpha
        #self.v_f_root = fan_2.v_f_root 
        #self.V_hor = fan_2.V_hor 
        self.V_c_slow = None 
        self.V_c_fast = None 

        self.related_fan1 = related_fan1
        self.related_fan2 = related_fan2

    def calc_V_c_slow(self):
        V_c_slow = (self.P_a / (self.k_v_f * self.related_fan2.mto_weight)) - ((self.V) * (self.related_fan2.rotor_alpha) + self.related_fan2.v_f_root)
        self.V_c_slow = V_c_slow
        return self.V_c_slow
  
    def calc_V_c_fast(self):
        V_c_fast = (self.P_a - (self.related_fan2.V_hor * self.D_h0 + self.related_fan2.T * self.related_fan2.v_f_root)) / (self.k_v_f * self.related_fan2.mto_weight)
        self.V_c_fast = V_c_fast
        return V_c_fast

fan_1 = Ducted_Fan_1(mass=float(718/4))
fan_2 = Ducted_Fan_2(mass=float(718/4), Cd0=0.05, V= 3, related_fan=fan_1)
fan_3 = Ducted_Fan_3(mass=float(718/4), gamma=0, related_fan1=fan_1, related_fan2 = fan_2 )

