import numpy as np
class Blade_Element_Theory:
    def __init__(self, theta_r, V_c, v_r, omega_r, c_r, d_r, rpm, r):
        #Gotta give the value when i create instance
        self.V_c = V_c
        self.v_r = v_r
        self.omega_r = omega_r
        self.theta_r = theta_r
        self.c_r = c_r
        self.d_r = d_r 
        self.rpm = rpm
        self.r = r 
    

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

    def calc_phi_1(self):
        self.phi_1 = np.arctan(self.V_c/self.omega_r)
        return self.phi_1
    
    def calc_phi_2(self):
        self.phi_2 = np.arctan(self.v_r/self.omega_r)
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
    
    