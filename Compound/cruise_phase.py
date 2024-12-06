# By Lintong
import numpy as np
import matplotlib.pyplot as plt 
# from ducted_fan_momentum_UI import func_min_locator

class WingedFlight:
    def __init__(self, vel, mass, chord, span, Cd0 = 0.022, power_a = None, wing_count=1, density = 1.225, g0 = 9.80665): #Torenbeek twin-engine piston 0.022 cd0
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

class LoiterFlight:
    def __init__(self):
        pass

class TransitionFlight:
    def __init__(self, mass, span, chord, Cl = 1.5, Cd0=0.01, TWR = 1, dt=0.001, gain_tilt =1, g0 = 9.80665): #assuming constant Cl
        self.g0 = g0
        self.mass = mass
        self.thrust_tot = self.mass * self.g0 * TWR
        self.L_win = 0
        self.D_win = 0
        self.dt = dt
        self.Cd0 = Cd0
        self.x = 0
        self.y = 30
        self.V_hor = 0.1
        self.V_ver = 0
        self.acc_hor = 0
        self.acc_ver = 0
        self.rotor_tilt = self.calc_RotorTilt(k=gain_tilt)

        
        self.Cl = Cl
        self.dyn_press = self.calc_DynPress()
        self.wing = Wing(dyn_press=self.dyn_press, span=span, weight=self.mass*self.g0, chord=chord, Cl=self.Cl)
        self.Cd = Cd0
        
        self.time = 0
        self.fd_x = []
        self.fd_y = []
        self.fd_V_hor = []
        self.fd_V_ver = []
        self.fd_acc_hor = []
        self.fd_acc_ver = []
        self.fd_rot_tilt = []
        self.fd_t = []
        self.fd_L_rot = []
        self.fd_L_win = []
        
        
        Run = True
        while Run:
            if self.L_win >= self.mass * self.g0:
                Run = False
                print('works')
            if self.y <= 0:
                Run = False
                print('crashed')
            #Updating Flight data
            self.time += dt 
            self.fd_x.append(self.x)
            self.fd_y.append(self.y)
            self.fd_V_hor.append(self.V_hor)
            self.fd_V_ver.append(self.V_ver)
            self.fd_acc_hor.append(self.acc_hor)
            self.fd_acc_ver.append(self.acc_ver)
            self.fd_rot_tilt.append(self.rotor_tilt)
            self.fd_t.append(self.time)
            #Finding Forces
            L_rot, T_rot = self.calc_RotorForces()
            self.L_win, self.D_win = self.calc_WingForces()
            self.fd_L_rot.append(L_rot)
            self.fd_L_win.append(self.L_win)
            # print(f'L_rot {L_rot}')
            # print(f'L_win {self.L_win}')
            # print(f'rotor_tilt {self.rotor_tilt}')
            #Updating Physics
            self.acc_hor, self.acc_ver = self.calc_Acceleration(lift=(L_rot+self.L_win), thrust=T_rot, drag=self.D_win)
            self.x, self.y, self.V_hor, self.V_ver = self.calc_UpdatePhysics(acc_hor=self.acc_hor, acc_ver=self.acc_ver)
            self.dyn_press = self.calc_DynPress()
            self.Cd = self.wing.calc_DragCoe(Cl=self.Cl)[0]
            # print(f'V_hor {self.V_hor}')
            # print(f'Cd {self.Cd}')
            #Updating AC Config
            self.rotor_tilt = self.calc_RotorTilt(k=gain_tilt)
            # print('-----------------')

        
    def calc_RotorForces(self): #keep hover thrust for entire phase, assuming constant fuel mass through phase as mass_ac >> fuelflow
        L_rotor = np.cos(self.rotor_tilt) * self.thrust_tot
        T_rotor = np.sin(self.rotor_tilt) * self.thrust_tot
        return L_rotor, T_rotor
    
    def calc_RotorTilt(self, k=1):
        rotor_tilt = k * self.V_hor
        if rotor_tilt > 90:
            rotor_tilt = 90
        rotor_tilt = rotor_tilt*(np.pi/180)
        return rotor_tilt
    
    def calc_Acceleration(self, lift, thrust, drag):
        force_forward_tot = thrust - drag
        acc_hor = force_forward_tot / self.mass
        force_upward_tot = lift - (self.mass * self.g0)
        acc_ver = force_upward_tot / self.mass
        return acc_hor, acc_ver
    
    def calc_UpdatePhysics(self, acc_hor, acc_ver):
        x_1 = self.x + self.V_hor * self.dt
        y_1 = self.y + self.V_ver * self.dt
        V_hor_1 = self.V_hor + acc_hor * self.dt
        V_ver_1 = self.V_ver + acc_ver * self.dt
        return x_1, y_1, V_hor_1, V_ver_1
    
    def calc_WingForces(self):
        L = self.Cl * self.dyn_press * self.wing.area
        D = self.Cd * self.dyn_press * self.wing.area
        return L, D
    
    def calc_DynPress(self):
        dyn_press = 0.5 * 1.225 * self.V_hor**2
        return dyn_press
    
class Wing:
    def __init__(self, dyn_press, span, weight, chord, Cl = 0, oswald_factor = 0.85, Cd0 = 0.008): #Torenbeek twin-engine piston 0.85 e 
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
        self.Cl = Cl
        self.Cd_tot = None
        self.Cd_ind = None
        
        Planform = self.calc_Planform()
        self.aspect_ratio = Planform[0]
        self.taper = Planform[1]
        self.area = Planform[2]
        
        if Cl == 0:
            self.Cl = self.calc_LiftCoe()
        
        drag_coe = self.calc_DragCoe(Cl=self.Cl)
        self.Cd_tot = drag_coe[0]
        self.Cd_ind = drag_coe[1]
        
        
    def calc_Planform(self):
        taper = self.chord_tp/self.chord_rt
        area = (self.chord_rt + self.chord_tp) * (self.span/2)
        aspect_ratio = self.span**2/area
        return aspect_ratio, taper, area
    
    def calc_LiftCoe(self):
        if self.dyn_press != 0:
            Cl = self.weight / (self.dyn_press*self.area)
            return Cl
    
    def calc_DragCoe(self, Cl): #Coefficient of drag
        Cd_ind = Cl**2 / (np.pi*self.aspect_ratio*self.oswald_factor)
        Cd_tot = self.Cd0 + Cd_ind
        return Cd_tot, Cd_ind
    
class Airfoil:
    def __init__(self, CL_max, Cd0,Cl0,Cl0_des,Cd0_des):
        self.CL_max=CL_max
        self.Cd0=Cd0
        self.Cl0=Cl0
        self.Cl0_des=Cl0_des
        self.Cd0_des=Cd0_des


def func_min_locator(list1, list2): #find the minimum and its index in list2 and locate the item at the same index in list1
    min_list2 = np.min(list2)
    min_list2_index = list(list2).index(min_list2)
    return list1[min_list2_index], min_list2

TEST = True
if TEST:
    #currently running B-29 root airfoil
    vel = np.arange(20,40,0.01)
    power_tot = []
    power_ind = []
    power_par = []
    Cl = []
    for i in vel:
        flight_point = WingedFlight(vel=i, power_a=10000, wing_count=2, mass= 1069.137, span=10, chord=[1]) #eq wing
        power_tot.append(flight_point.power_tot)
        power_ind.append(flight_point.power_ind)
        power_par.append(flight_point.power_par)
        Cl.append(flight_point.wing.Cl)
        
    power_tot = np.array(power_tot)/1000
    power_ind = np.array(power_ind)/1000
    power_par = np.array(power_par)/1000
    
    # print(func_min_locator(vel, power_tot))
    
    plt.plot(vel, power_tot, "-", label="Total")
    plt.plot(vel, power_ind, "--", label="induced")
    plt.plot(vel, power_par, "--", label="parasitic")
    
    vel_power_min, power_min = func_min_locator(vel, power_tot)
    plt.plot(vel_power_min, power_min, '.', label=f'Min. Power \nV = {vel_power_min:.1f}m/s \nP = {power_min:.2f}kW')
    
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Power Required [kW]')
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.clf
    
    plt.plot(vel, np.array(Cl), "-", label="Cl")
    plt.axhspan(1.5, max(np.array(Cl)), color='red', alpha=0.5)
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('C_l Required [-]')
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.clf
    
    #Transition Flight
    tr = TransitionFlight(mass=718, span=10, Cl=1.5, chord=[1], TWR=1.0, gain_tilt=3.2)
    #gain = 3 for minimum alt change
    plt.plot(tr.fd_t, tr.fd_L_rot, "-", label = 'Lift from rotor')
    plt.plot(tr.fd_t, tr.fd_L_win, "-", label = 'Lift from wing')
    plt.plot(tr.fd_t, np.array(tr.fd_L_win)+np.array(tr.fd_L_rot), "-", label = 'Total Lift')
    # plt.title('')
    plt.legend()
    plt.grid(True)
    plt.xlabel('Time [s]')
    plt.ylabel('Lift [s]')

    # plt.plot(tr.fd_t, tr.fd_rot_tilt, "-")
    # plt.title('Flight Trtrajectory')
    # plt.xlabel('Downrange [m]')
    # plt.ylabel('Height [m]')
    
    plt.show()
    # plt.clf