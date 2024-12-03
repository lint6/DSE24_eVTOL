# By Lintong
import numpy as np
import matplotlib.pyplot as plt 

class WingedFlight:
    def __init__(self, vel, mass, chord, span, Cd0 = 0.1, power_a = None, wing_count=1, density = 1.225, g0 = 9.80665):
        self.power_a = power_a #power available
        self.wing_count = wing_count    #amount of wing on aircraft
        self.vel = vel  #airspeed of aircraft
        self.mass = mass    #mass of aircraft
        self.weight = self.mass*g0  #weight of aircraft
        self.density = density  #density of env
        self.Cd0_nw = Cd0   #list of Cd0 of non-wing parts
        self.chord = chord  #chord of wings
        self.span = span    #span of wings
        
        self.weight_div = self.calc_MassDistrib() #Distribution of lift generation between wings
        self.dyn_press = self.calc_DynPress()
        self.wing = Wing(dyn_press=self.dyn_press, span=self.span, weight=self.weight_div, chord=self.chord) #Create the wing object
        
        Cd_calc = self.calc_DragCoe()
        self.Cd = Cd_calc[0]
        self.Cd0 = Cd_calc[1]
        self.Cd_ind = Cd_calc[2]
        
        drag = self.calc_DragForce()
        self.drag_tot = drag[0]
        self.drag_ind = drag[1]
        self.drag_par = drag[2]
        
        power = self.calc_PowerRequiredLevel()
        self.power_tot = power[0]
        self.power_ind = power[1]
        self.power_par = power[2]
        
        
    def calc_PowerRequiredLevel(self):
        power_tot = self.vel * self.drag_tot
        power_ind = self.vel * self.drag_ind
        power_par = self.vel * self.drag_par
        return power_tot, power_ind, power_par
    
    def calc_DragCoe(self): #Coefficient of drag
        Cd0_sum = np.sum(self.Cd0_nw)
        Cd_ind_sum = 0
        Cd0_sum += self.wing.Cd0 * self.wing_count
        Cd_ind_sum += self.wing.Cd_ind * self.wing_count
        Cd_tot = Cd0_sum + Cd_ind_sum
        return Cd_tot, Cd0_sum, Cd_ind_sum
    
    def calc_DragForce(self):
        drag_force = self.Cd * self.dyn_press * self.wing.area
        drag_force_ind = self.Cd_ind * self.dyn_press * self.wing.area
        drag_force_par = self.Cd0 * self.dyn_press * self.wing.area
        return drag_force, drag_force_ind, drag_force_par
    
    def calc_ClimbRate(self):
        pass
    
    def calc_PowerSurplus(self): #Difference in power required and power availiable
        pass
    
    def calc_DynPress(self):
        dyn_press = 0.5 * self.density * self.vel**2
        return dyn_press
    
    def calc_MassDistrib(self):
        return self.weight/self.wing_count
    
class Wing:
    def __init__(self, dyn_press, span, weight, chord, oswald_factor = 0.81, Cd0 = 0.008):
        self.span = span

        self.weight = weight
        self.chord = chord #list of chord, single input means constant chord
        
        self.dyn_press = dyn_press
        self.Cd0 = Cd0
        self.oswald_factor = oswald_factor
        
        self.chord_rt = chord[0]
        self.chord_tp = chord[-1]
        
        self.aspect_ratio = None
        self.taper = None
        self.area = None
        self.Cl = None
        self.Cd_tot = None
        self.Cd_ind = None
        
        Planform = self.calc_Planform()
        self.aspect_ratio = Planform[0]
        self.taper = Planform[1]
        self.area = Planform[2]
        
        self.Cl = self.calc_LiftCoe()
        
        drag_coe = self.calc_DragCoe()
        self.Cd_tot = drag_coe[0]
        self.Cd_ind = drag_coe[1]
        
        
    def calc_Planform(self):
        taper = self.chord_tp/self.chord_rt
        area = (self.chord_rt + self.chord_tp) * (self.span/2)
        aspect_ratio = self.span**2/area
        return aspect_ratio, taper, area
    
    def calc_LiftCoe(self):
        Cl = self.weight / (self.dyn_press*self.area)
        return Cl
    
    def calc_DragCoe(self): #Coefficient of drag
        Cd_ind = self.Cl**2 / (np.pi*self.aspect_ratio*self.oswald_factor)
        Cd_tot = self.Cd0 + Cd_ind
        return Cd_tot, Cd_ind
    
class Airfoil:
    def __init__(self, CL_max, Cd0):
        self.CL_max=CL_max
        self.Cd0=Cd0
        
TEST = True
if TEST:
    vel = np.arange(10,60,0.1)
    power_tot = []
    power_ind = []
    power_par = []
    Cl = []
    for i in vel:
        flight_point = WingedFlight(vel=i, power_a=10000, wing_count=2, mass= 718, span=10, chord=[1.25])
        power_tot.append(flight_point.power_tot)
        power_ind.append(flight_point.power_ind)
        power_par.append(flight_point.power_par)
        Cl.append(flight_point.wing.Cl)
        
    power_tot = np.array(power_tot)/1000
    power_ind = np.array(power_ind)/1000
    power_par = np.array(power_par)/1000
    plt.plot(vel, power_tot, "-", label="Total")
    plt.plot(vel, power_ind, "-", label="induced")
    plt.plot(vel, power_par, "-", label="parasitic")
    plt.legend()               
    plt.grid(True)
    plt.show()
    plt.clf
    
    plt.plot(vel, np.array(Cl), "-", label="Cl")
    plt.grid(True)
    plt.show()