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
        Cd0_sum += self.wing.Cd0
        Cd_ind_sum += self.wing.Cd_ind
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

class SteadyTransitionFlight:
    def __init__(self, V, mass, span, chord, radius, wing_count = 1, Cd0=0.022, Cl = 1.5):
        self.mass = mass
        self.g0 = 9.80665
        self.weight = self.mass * self.g0
        self.V = V
        self.wing = Wing(dyn_press=0, span=span, weight=self.weight, chord=chord, Cl=Cl)
        self.Cd0_fus = Cd0
        self.radius = radius
        self.wing_count = wing_count
        self.dyn_press = self.calc_DynPress()
        
        self.WingForce = self.calc_WingForces() #0 is lift, 1 is drag
        
        self.T_hor = self.WingForce[1]
        self.T_ver = self.weight - self.WingForce[0]
        
        self.rotor = self.calc_RotorCon() #0 is angle, 1 is thrust
        self.power = self.calc_Power() #0 is total, 1 is vertical, 2 is horizontal
        
    def calc_DynPress(self):
        dyn_press = 0.5 * 1.225 * self.V**2
        return dyn_press
    
    def calc_WingForces(self): #from transition flight code
        Lift = self.wing.Cl * self.dyn_press * self.wing.area * self.wing_count
        Drag = (self.wing.Cd_tot + self.Cd0_fus) * self.dyn_press * self.wing.area
        return Lift, Drag
    
    def calc_RotorCon(self): #Rotor condition, angle of rotor and thrust produced
        if self.T_hor:
            angle = np.arctan(self.T_ver / self.T_hor)
        else:
            angle = np.pi/2
        thrust = (self.T_ver**2 + self.T_hor**2)**(1/2)
        return angle, thrust
    
    def calc_Power(self, density = 1.225):
        power_ver = self.T_ver * np.sqrt(self.T_ver/(2* np.pi * self.radius**2 * density*4))
        power_hor = self.T_hor * self.V
        power_tot = power_hor + power_ver
        return power_tot, power_ver, power_hor
        
    
    
class TransitionFlight:
    def __init__(self, mass, span, chord, Cl = 1.5, Cd0=0.01, TWR = 1, dt=0.001, gain_tilt =1, g0 = 9.80665, Run = True): #assuming constant Cl
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
        
        Run = Run
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


#This is the graph to add the tnagent line
#here
#here
#here
#here 
TEST = True
if TEST:
    chord = [1]
    span = 10
    mass = 1069.137
    rotor_r = 0.9652/2
    wing_num = 2 #doubles the lift, doesnt use equivalent wing 
    # currently running B-29 root airfoil
    vel = np.arange(20,60,0.001)
    power_tot = []
    power_ind = []
    power_par = []
    Cl = []
    for i in vel:
        flight_point = WingedFlight(vel=i, power_a=10000, wing_count=2, mass= 1069.137, span=10, chord=[1])
        power_tot.append(flight_point.power_tot)
        power_ind.append(flight_point.power_ind)
        power_par.append(flight_point.power_par)
        Cl.append(flight_point.wing.Cl)
        
    power_tot = np.array(power_tot)/1000
    power_ind = np.array(power_ind)/1000
    power_par = np.array(power_par)/1000
    
    tolerance = 0.001
    tangent = np.gradient(power_tot,vel)
    intersection_point = None

    for x, y, slope in zip(vel, power_tot, tangent):
        if x != 0:
            slope_from_origin = y/x
            if abs(slope-slope_from_origin) <= tolerance:
                intersection_point = (x,y)
                break

    # print(func_min_locator(vel, power_tot))
    
    plt.plot(vel, power_tot, label="Total Power", color = "red", linewidth = 2)
    plt.plot(vel, power_ind, label="Induced Power", color = "green")
    plt.plot(vel, power_par, label="Parasitic Power", color = "black")
    
    vel_power_min, power_min = func_min_locator(vel, power_tot)
    plt.plot(vel_power_min, power_min, '.', label=f'Min. Power \nV = {vel_power_min:.1f}m/s \nP = {power_min:.2f}kW')

    print(intersection_point)

    if intersection_point:
        x_tangent, y_tangent = intersection_point
        tangent_slope = y_tangent / x_tangent
        tangent_x = np.linspace(0, max(vel), 100)
        tangent_y = tangent_slope * tangent_x
        plt.scatter([x_tangent], [y_tangent], color = "orange", label = "max range")
        plt.plot(tangent_x, tangent_y, linestyle = "dotted", linewidth = 1, color = "black")
    
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
    plt.margins(0)
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.clf
    
    # #Transition Flight
    # tr = TransitionFlight(mass=mass, span=span, Cl=1.5, chord=chord, TWR=1.0, gain_tilt=2.7)
    # #gain = 3 for minimum alt change
    # plt.plot(tr.fd_t, np.array(tr.fd_L_rot)/1000, "-", label = 'Lift from rotor')
    # plt.plot(tr.fd_t, np.array(tr.fd_L_win)/1000, "-", label = 'Lift from wing')
    # plt.plot(tr.fd_t, (np.array(tr.fd_L_win)+np.array(tr.fd_L_rot))/1000, "-", label = 'Total Lift')
    # # plt.title('')
    # plt.legend()
    # plt.grid(True)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Lift [kN]')
    # plt.show()
    # plt.clf

    # plt.plot(tr.fd_t, tr.fd_rot_tilt, "-")
    # # plt.title('Flight Trtrajectory')
    # plt.xlabel('Downrange [m]')
    # plt.ylabel('Height [m]')
    
    # plt.show()
    # plt.clf
    
    '''Transition Flight'''
    vel1 = np.arange(0,60,0.01)
    thrust = []
    angle = []
    thrust_ver = []
    thrust_hor = []
    power_tr = []
    power_tr_v = []
    power_tr_h = []
    for i in vel1:
        plane = SteadyTransitionFlight(V=i, mass=mass, span=span, chord=chord, radius=rotor_r, wing_count=wing_num, Cl = 1.5)
        if plane.power[1] > 0:
            thrust.append(plane.rotor[1])
            angle.append(plane.rotor[0])
            thrust_ver.append(plane.T_ver)
            thrust_hor.append(plane.T_hor)
            power_tr.append(plane.power[0])
            power_tr_v.append(plane.power[1])
            power_tr_h.append(plane.power[2])
        
    vel1 = list(vel1)[:len(thrust)]
    
    '''Winged flight'''
    vel2 = np.arange(vel1[-1],130,0.01)
    power_tot = []
    power_ind = []
    power_par = []
    Cl = []
    for i in vel2:
        flight_point = WingedFlight(vel=i, power_a=10000, wing_count=wing_num, mass= mass, span=span, chord=chord) #eq wing
        power_tot.append(flight_point.power_tot)
        power_ind.append(flight_point.power_ind)
        power_par.append(flight_point.power_par)
        Cl.append(flight_point.wing.Cl)


    plt.plot(vel1, np.array(power_tr)/1000, "-", label = "Transition Power")
    plt.plot(vel1, np.array(power_tr_h)/1000, "--", label = "Transition Horizontal Power")
    plt.plot(vel1, np.array(power_tr_v)/1000, "--", label = "Transition Vertical Power")
    plt.plot(vel2, np.array(power_tot)/1000, "-", label = "Total Power")
    plt.plot(vel2, np.array(power_ind)/1000, "--", label = "Induced power")
    plt.plot(vel2, np.array(power_par)/1000, "--", label = "Parasite Power")
    

    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Power Required [kW]')
    plt.legend()
    plt.grid(True)
    plt.margins(0.01)
    plt.show()
    plt.clf
    
    plt.plot(vel1, np.array(angle)/np.pi*180)
    plt.show()

    """
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
    """

    """
    tolerance = 1e-2
    tangent = np.gradient(power_tot,vel1)
    intersection_point = None
    for x, y, slope in zip(vel1, power_tot, tangent):
        if x != 0:
            slope_from_origin = y/x
            if abs(slope-slope_from_origin) <= tolerance:
                intersection_point = (x,y)
                break
"""