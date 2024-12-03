import numpy as np
from co_rotor_sizing import RotorSizing

class SoundAnalysis:
    def __init__(self):
        #Conversions
        self.m_to_f = 1/0.3048 #Conversion factor from meter to feet
        self.N_to_lbs = (1/9.80665)*2.20462

        #Input
        self.rotorsizing = RotorSizing()

        #Observer position relative to rotor
        self.x = 0 #np.linspace(0,100,1000) #Measured in direction of motion of helicopter [ft]
        self.y = 0 #np.linspace(-100,100,1000) #Measured at 90 deg to x in disc plane [ft]
        self.z = 500 #Flyover height [ft]
        self.r = (self.x**2+self.y**2+self.z**2)**0.5

        #Inputs for rotational noise
        self.R = self.rotorsizing.rotor_radius*self.m_to_f #Get rotor radius in [f]
        self.A = np.pi*(self.R**2) #Rotor area [ft^2]
        self.n = self.rotorsizing.omega #Rotor rotational speed [rad/s] 
        self.V = self.rotorsizing.max_forward_velocity*self.m_to_f #Set to max speed for now [ft/s]
        self.c = self.rotorsizing.speed_of_sound*self.m_to_f #Speed of sound [ft/s]
        self.B = self.rotorsizing.n_blades
        self.T = self.rotorsizing.T_forward_flight*self.N_to_lbs #Thrust [lbs]

        #Inputs for vortex noise
        self.D = 2*self.R #Rotor diameter [ft]
        self.V_07 = 0.7*((self.n*np.pi*self.D)/60) #Linear speed of 0.7 radius section
        self.chord = self.rotorsizing.chord*self.m_to_f #Chord [ft]
        self.A_b = self.B*self.R*self.chord #Total blade area [ft^2]

        #Initialize
        self.rotational_noise()
        self.vortex_noise()

    def rotational_noise(self):
        #Calculate rotational Mach number
        self.M = (0.8*self.n*self.R)/(self.c)

        #Calculate flight Mach number
        self.M_f = self.V/self.c

        #Calculate theta dash
        self.theta_dash = np.arccos(self.x/self.r)

        #Calculate effective rotational Mach number
        self.M_E = self.M/(1-self.M_f*np.cos(self.theta_dash))

        #Calculate theta
        self.theta = np.pi/2 #Set to angle that gives max noise for now (unused)

        #Interpolate from graph of first harmonic
        M_series = np.linspace(0.2,1,9)
        noise_series = np.array([77,82,85,88,91,93,95,96,97])
        interpolator = np.polyfit(M_series, noise_series,3)

        self.rotational_SPL_uncorrected = interpolator[0]*(self.M_E)**3  + interpolator[1]*(self.M_E)**2 + interpolator[2]*(self.M_E) + interpolator [3]

        #Apply SPL correction
        self.rotational_SPL = self.rotational_SPL_uncorrected + 11 + 10*np.log10((self.T/(self.r**2))*(self.T/self.A))

        #Calculate fundamental frequency
        self.f_rotational = (self.n*self.B)/(2*np.pi*(1-self.M_f*np.cos(self.theta)))
    def vortex_noise(self):
        #Calculate SPL for 300 ft distance
        self.vortex_SPL_uncorrected = 10*(2*np.log10(self.V_07)+2*np.log10(self.T)-np.log10(self.A_b)-3.57)

        #Correct for distance (500 ft)
        self.vortex_SPL = self.vortex_SPL_uncorrected - 20*np.log10(self.z/300)

        #Calculate frequency
        self.thickness = 0.12*self.chord #NACA0012 airfoil thickness to chord
        self.f_vortex = (self.V_07*0.28)/self.thickness  #Only valid for small AoA

    def display_parameters_rotor(self):
        print(f"Equivalent M: {self.M_E}")
        print(f"Uncorrected SPL: {self.rotational_SPL_uncorrected} dB")
        print(f"Correction factor: {self.rotational_SPL-self.rotational_SPL_uncorrected} dB")
        print(f"Corrected SPL: {self.rotational_SPL} dB")
        print(f"Fundamental frequency: {self.f_rotational} Hz")
    def display_paramenters_vortex(self):
        print(f"Uncorrected SPL: {self.vortex_SPL_uncorrected} dB")
        print(f"Corrected SPL: {self.vortex_SPL} dB")
        print(f"Vortex frequency: {self.f_vortex} Hz")


if __name__ == '__main__':
    sound = SoundAnalysis()
    print('-----Rotational noise-----')
    sound.display_parameters_rotor()
    print('-----Vortex noise-----')
    sound.display_paramenters_vortex()

