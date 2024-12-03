import numpy as np 
from math import *
import scipy 
import matplotlib.pyplot as plt 
import sys
#Code from YUNJAE 


#def calc_p_idd_kappa(self):
#    self.calc_power_ideal_hover()
#    self.calc_kappa_d()
#    self.p_idd_kappa = self.P_id_h * self.kappa_d
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
    def __init__(self, mass, radius=None, TWR= 1 , V_c=None, density=1.225, g0 = 9.80665, P_a = None): 
        # if V_c== None and P_a == None:
        #     raise Exception('Vertical rate and power are both none, need at least one')
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.g0 = g0

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

    # def calc_thrust_loading(self):
    #     self.calc_disc_loading()
    #     self.thrust_loading = self.disc_loading * self.TWR * 9.81
    #     return self.thrust_loading

    def calc_v_h(self): #Induced velocity during hover
        self.calc_disc_loading()
        self.v_h = np.sqrt(self.thrust_loading/(2 * self.density)) #eq 2.13
        return self.v_h

    def calc_thrust_hover(self): #Thrust required for hover
        self.calc_v_h()
        self.thrust = 2 * (np.pi * self.radius**2 ) * self.density * (self.v_h)**2 #eq 2.12
        return self.thrust
    
    def calc_power_ideal_hover(self): #Power require for hover
        self.calc_mass()
        self.calc_v_h()
        self.P_id_h = self.weight * self.v_h
        return self.P_id_h

    def calc_v_c (self): #Induced velocity during axial maneuver
        self.calc_v_h()
        if self.V_c == None:
            self.calc_V_c_kappa()
        else:
            self.calc_V_c_nd()
        self.v_c_nd = -0.5 * self.V_c_nd + np.sqrt(0.25 * self.V_c_nd**2 + 1) #eq 2.14a
        self.v_c = self.v_c_nd * self.v_h
        return self.v_c, self.v_c_nd
    
    def calc_kappa_p(self): #kappa from known power, unknown climb rate
        self.calc_power_ideal_hover()
        self.kappa = self.P_a / self.P_id_h
        return self.kappa
    
    def calc_kappa_Vc(self): #kappa from known climb rate, unknown power
        self.calc_power_ideal_hover()
        self.calc_V_c_nd()
        self.kappa = 0.5 * self.V_c_nd + sqrt(0.25 * self.V_c_nd ** 2 + 1)
        return self.kappa
    
    def calc_V_c_kappa(self): #Axial rate of climb given kappa
        self.calc_kappa_p()
        self.calc_v_h()
        self.V_c_nd = self.kappa - (1/self.kappa)
        self.V_c = self.V_c_nd * self.v_h
        return self.V_c, self.V_c_nd
    
    def calc_V_c_nd(self): #Find V_c_nd if calc_V_c_kappa() is not used
        self.calc_v_h()
        self.V_c_nd = self.V_c / self.v_h
        return self.V_c_nd
    
    def calc_P_kappa(self): #Axial power required given kappa
        self.calc_kappa_Vc()
        self.calc_v_h()
        self.calc_mass()
        self.P_id_c = self.weight * self.v_h * self.kappa
        return self.P_id_c

    '''
    Commented due to Vd/vh>2 is unlikely to encounter during normal operation
    # def calc_kappa_d(self):
    #     self.calc_V_c_kappa()
    #     self.calc_v_d()
    #     if np.abs(self.V_c_kappa/self.v_d)>=2 :
    #         self.kappa_d = 0.5 * (self.V_c_kappa) + np.sqrt(0.25 * self.V_c_kappa**2 - 1)
    #     elif np.abs(self.V_c_kappa/self.v_d)<1 :
    #         self.kappa_d = 0.5 * (self.V_c_kappa) + np.sqrt(0.25 * self.V_c_kappa**2 + 1)
    #     return self.kappa_d
    #
    # def calc_p_idd(self):
    #     self.calc_V_c_kappa()
    #     self.calc_v_d()
    #     if np.abs(self.V_c_kappa/self.v_d)>=2 :
    #         self.p_idd = - self.mass * 9.81 * (self.v_d + self.V_c_kappa)
    #     elif np.abs(self.V_c_kappa/self.v_d)<1 :
    #         self.p_idd = self.mass * 9.81 * (self.v_d + self.V_c_kappa)
    #     return self.p_idd
    #
    # def calc_V_d_kappa_d(self):
    #     self.calc_V_c_kappa()
    #     self.calc_v_d()
    #     if np.abs(self.V_c_kappa/self.v_d)>=2 :
    #         V_d_kappa_nd = - (self.kappa_d + (1/self.kappa_d))
    #         self.V_d_kappa_d = V_d_kappa_nd * self.v_h  
    #     elif np.abs(self.V_c_kappa/self.v_d)<1 :
    #         V_d_kappa_nd = - self.kappa_d + (1/self.kappa_d)
    #         self.V_d_kappa_d = V_d_kappa_nd * self.v_h
    #     return self.V_d_kappa_d
    '''

    
class Ducted_Fan_2: #pure horizontal
    def __init__(self, mass, Cd0, V, radius=0.625, TWR= 1 , gamma = 0, V_c=None, density=1.225, k_v_f = 1, related_fan = None): 
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.Cd0 = Cd0
        self.k_v_f = k_v_f 
        self.V = V # we need to know this somehow first .... 
        self.gamma = gamma

        # Optional Inputs with Defaults
        self.TWR = TWR  # Thrust-to-weight ratio
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
    def __init__(self, mass, gamma, radius=0.625, TWR= 1 , V_c=None, V =4, density=1.225, D_h0 = 500, k_v_f = 1, P_a =45000,  related_fan1 = None, related_fan2 = None): 
        # Required Inputs
        self.mass = mass  # Maximum takeoff weight
        self.radius = radius  # Fan radius (m)
        self.D_h0 = D_h0 #Drag horizontal
        self.k_v_f = k_v_f  #coefficient of vertical drag ~ 1
        self.V = V # Velocity we give
        self.gamma = gamma # airpath angle
        self.P_a = P_a # power avilable

        # Optional Inputs with Defaults
        self.TWR = TWR  # Thrust-to-weight ratio
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


