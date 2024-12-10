import numpy as np 
from math import *
import scipy 
import matplotlib.pyplot as plt 
import sys
# cl limit +12 to - 15 
# airfoil cd0 is 0.008 

def find_roots(coefficients):
    if len(coefficients) != 5:
        raise ValueError("You must provide exactly 5 coefficients for a 4th order polynomial.")
    roots = np.roots(coefficients)  
    real_roots = roots[np.isreal(roots)]
    real_roott = min(real_roots[real_roots>0])
    real_root = real_roott.real
    return real_root


class Ducted_Fan: #angled flight
    def __init__(self, mass, gamma, V, Cd0, k_v_f, radius, TWR= 1 , V_c=None, density=1.225, g0 = 9.80665, P_a = None): 
        # if V_c== None and P_a == None:
        #     raise Exception('Vertical rate and power are both none, need at least one')
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.g0 = g0
        self.Cd0 = Cd0 # 0.022 assumed airplane fuselage drag coeff
        self.gamma = gamma 
        self.V = V 
        self.k_v_f = k_v_f


        # Optional Inputs with Default
        self.TWR = TWR  # Thrust-to-weight ratio
        self.V_c = V_c  # Climb freestream velocity (m/s)
        self.density = density  # Air density (kg/m^3)
        self.P_a = P_a # Power available set by user

        # Calculated Attributes (initialized as None)
        self.V_c_nd = None
        self.weight = None
        self.disc_loading = None
        self.thrust_loading = None 
        self.v_h = None  # Hover induced velocity
        self.thrust = None # thurst
        #self.hover_induced_velocity = None
        self.v_d = None #induced velocitt decent
        self.v_d_nd = None #above non dimensional
        self.v_c = None #induced velocity climb
        self.v_c_nd = None #above non dimensional
        self.kappa = None # #ratio of ideal power abaliable o ideal power required in hover, from climb rate
        self.P_id_c = None # power ideal climb
        self.P_id_h = None  # power ideal hover
        # self.kappa = None  # #ratio of ideal power abaliable o ideal power required in hover, from power avalibale 
        self.p_idd = None # power ideal decent
        self.kappa_d = None  # #ratio of ideal power abaliable o ideal power required in hover, from power avalable
        self.p_idd_kappa = None # power ideal descent from, kappa 
        self.V_d_kappa_d = None  #vertical Descnet rate from kappa RATE


    def calc_mass(self):
        new_mass = -sys.maxsize
        if new_mass > self.mass:
            self.mass = new_mass
        self.weight = self.mass * self.g0
        return self.mass, self.weight
    
    def calc_disc_loading(self):
        self.calc_mass()
        self.disc_loading = self.mass / (np.pi * self.radius**2 ) #Definition of disc loading
        self.thrust_loading = self.weight * self.TWR / (np.pi * self.radius**2 ) #Definition of thrust area loading
        return self.disc_loading, self.thrust_loading

    def calc_v_h(self): #Induced velocity during hover
        self.calc_disc_loading()
        self.v_h = np.sqrt(self.thrust_loading/(2 * self.density)) #eq 2.13
        return self.v_h
    
    
    def calc_V_horizontal(self):
        self.V_hor = self.V * np.cos(self.gamma)
        return self.V_hor
    
    def calc_V_nd(self):
        self.calc_v_h()
        self.V_nd = self.V / self.v_h
        return self.V_nd
    
    def calc_D(self):
        self.D = (0.5) * self.density * (self.V **2) * self.Cd0 * (self.radius**2 * np.pi) # body parasite drag
        # iunduced drag
        #paraise drage

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
        self.v_f_root = find_roots(coefficients=my_coefficient) * self.v_h
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
        return self.p_idf #power_induced, power_parasitic

#     def __init__(self, mass, gamma, V, Cd0, k_v_f, radius=None, TWR= 1 , V_c=None, density=1.225, g0 = 9.80665, P_a = None): 

fan_exp = Ducted_Fan(mass=1069.137, gamma= 0, V= 50, Cd0=1, k_v_f=1, radius=(0.9652/2))
print(fan_exp.calc_p_idf())
