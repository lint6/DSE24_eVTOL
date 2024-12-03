import numpy as np
# Code is made of YUNJAE 
class Blade_Element_Theory:
    def __init__(self, theta_r, theta_0, theta_tot,  V_c, v_r, omega, sigma, c_r, d_r, rpm, r, R, T, Q, P):
        #Gotta give the value when i create instance
        self.V_c = V_c
        self.v_r = v_r
        self.omega = omega
        self.sigma = sigma 
        self.c_r = c_r
        self.d_r = d_r 
        self.rpm = rpm
        self.r = r 
        self.R = R
        self.theta_r = theta_r
        self.theta_0 = theta_0
        self.theta_tot = theta_tot 

        #Check where this comes from, maybe we need to redefine the T, Q, P as we are giving it as input vairbale here 
        self.T = T 
        self.Q = Q 
        self.P = P 

    


        #fixed input variable which appplies to all the class def and no iterating
        self.a = 10 
        self.density = 1.225
        self.Cd_r = 1


        self.alpha_blade = None 
        self.phi_1 = None 
        self.phi_2 = None 
        self.phi = None 
        self.U_r = None 
        self.dL_r = None 
        self.dD_p_rr = None 
        self.dD_p_rr_simple = None 
        self.dT_r = None 
        self.dQ_r = None 
        self.dP_r = None 
        self.dT_r1 = None 
        self.theta_r_bar = None 
        self.V_t = None 
        self.v_r_bar_nd = None 
        self.alpha_r_bar = None 
        self.c_l_r_bar = None 

        self.C_T = None 
        self.C_Q = None 
        self.C_P = None 

        self.c_t = None 
        self.c_q = None 
        self.c_p = None 

    def calc_phi_1(self):
        self.phi_1 = np.arctan(self.V_c/(self.omega * self.r))
        return self.phi_1
    
    def calc_phi_2(self):
        self.phi_2 = np.arctan(self.v_r/(self.omega * self.r))
        return self.phi_2
    
    def calc_phi(self):
        self.calc_phi_1()
        self.calc_phi_2()
        self.phi = self.phi_1 + self.phi_2
        return self.phi

    def calc_alpha_blade(self):
        self.alpha_blade = self.theta_r - (self.phi_1 + self.phi_2)
        return self.alpha_blade
    
    def calc_cl_r(self):
        self.cl_r = self.alpha_blade * self.a 
        return self.cl_r
    
    def calc_U_r(self):
        self.U_r = self.V_c + self.v_r + self.rpm * self.r
        return self.U_r

    def calc_dL_r(self):
        self.dL_r = 0.5 * self.density * self.a * self.alpha_blade*(self.U_r**2) * self.c_r * self.d_r
        return self.dL_r
    
    def calc_dD_p_rr(self):
        self.dD_p_rr = 0.5 * self.Cd_r * self.density * self.U_r * self.c_r * self.d_r
        return self.dD_p_rr
    
    def calc_dD_p_rr_simple(self):
        self.dD_p_rr_simple = 0.5 * self.Cd_r * self.density * ((self.rpm * self.r)**2) * self.c_r * self.d_r
        return self.dD_p_rr_simple
    
    def calc_dT_r(self):
        self.dT_r = self.dL_r * np.cos(self.phi) - self.dD_p_rr * np.sin(self.phi)
        return self.dT_r
    
    def calc_dQ_r(self):
        self.dQ_r = (self.dL_r * np.sin(self.phi) + self.dD_p_rr * np.cos(self.phi)) * self.r
        return self.dQ_r
    
    def calc_dP_r(self):
        self.dP_r = self.dQ_r * self.rpm
        return self.dP_r
    
    # Combined blade element theory with momentum theory 

    def calc_dT_r1(self):
        self.dT_r1 = 4 * np.pi * self.density * (self.V_c + self.v_r) * self.v_r * self.r * self.d_r
        return self.dT_r1
    
    # Twist function !!!!!!!!!!!!!!!!!
    def calc_theta_r_bar(self):
        self.theta_r_bar = self.theta_0 - self.theta_tot * self.r/self.R 
        return self.theta_r_bar
    
    def calc_V_t(self):
        self.V_t = self.omega * self.R 
        return self.V_t
    
    def calc_v_r_bar(self):
        self.v_r_bar = self.V_t * ((- (self.a * self.sigma / 16)) + np.sqrt((self.a * self.sigma / 16)**2 + self.a * self.sigma * (self.r/self.R)* self.theta_r_bar)/8)
        return self.v_r_bar
    
    def calc_v_r_bar_nd(self):
        self.v_r_bar_nd = self.v_r_bar / self.V_t
        return self.v_r_bar_nd

    #Angle of attack at a station r_bar 
    def calc_alpha_r_bar(self):
        self.alpha_r_bar = self.theta_0 - np.arctan(self.v_r_bar_nd/self.v_r_bar)
        return self.alpha_r_bar
    
    def calc_c_l_r_bar(self):
        self.c_l_r_bar = self.alpha_r_bar * self.a
        return self.c_l_r_bar
    

    #Non dimentional coefficient 
    def calc_thrust_coefficient(self):
        self.C_T = self.T / (np.pi * (self.R **2) * self.density * (self.V_t ** 2 ))
        return self.C_T
    
    def calc_torque_coefficient(self):
        self.C_Q = self.Q / (np.pi * (self.R **3) * self.density * (self.V_t ** 2 ))
        return self.C_Q
    
    def calc_power_coefficient(self):
        self.C_P = self.P / (np.pi * (self.R **2) * self.density * (self.V_t ** 3 ))
        return self.C_P
    
    # Non dimentional against the total blade area 

    def calc_local_thurst_coefficient(self):
        self.c_t = self.C_T / self.sigma
        return self.c_t
    
    def calc_local_torque_coefficient(self):
        self.c_q = self.C_Q / self.sigma
        return self.c_q

    def calc_local_power_coefficient(self):
        self.c_p = self.C_P / self.sigma
        return self.c_p
    
    