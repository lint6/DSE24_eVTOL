# Created by Lintong
import numpy as np
import matplotlib.pyplot as plt 

def calc_AircraftForce(V_free, Cl_max = 1.4, Cd0=0.022+0.008, area=10, AR=10, e = 0.85): #Torenbeek
    dyn_press = 1/2 * 1.225 * V_free**2 #only critcal condition
    L = Cl_max*dyn_press*area*2
    D = (Cd0 + Cl_max**2/(np.pi*e*AR))*dyn_press*area
    return L, D
def calc_MaxThrust(mass):
    weight = mass * 9.80665
    thrust = weight * 1.01
    return thrust, weight
def calc_AdvanceRatio(V_free, eta, radius=0.9652/2, RPM=600):
    V_tip = radius * RPM
    V_locl = V_free * np.cos(eta-np.pi/2)
    AR = V_locl / V_tip
    return AR
def calc_ThrustVector(thrust, eta):
    thrust_hor = thrust * np.cos(eta)
    thrust_ver = thrust * np.sin(eta)
    return thrust_hor, thrust_ver
def func_Envelope(mass, eta, V, adv_ratio_lim = 0.4):
    eta = np.radians(eta)
    adv_ratio = calc_AdvanceRatio(V_free=V, eta=eta)
    thrust, weight = calc_MaxThrust(mass=mass)
    thrust_x, thrust_y = calc_ThrustVector(thrust=thrust, eta=eta)
    lift, drag = calc_AircraftForce(V_free=V)
    excess_lift = lift - weight
    print(excess_lift)
    drag_con = thrust_x - drag

    
    adv_ratio = adv_ratio < adv_ratio_lim
    drag_con = drag_con >= 0
    lift_con = thrust_y + excess_lift >= 0
    grid_out = np.logical_and(adv_ratio, drag_con)
    grid_out = np.logical_and(grid_out, lift_con)
    return grid_out, drag_con, lift_con, adv_ratio, excess_lift, thrust_y
    
    
    
    
Run = True
if Run:
    mass = 1069.137
    V_max = 150
    V_step = 0.1
    eta_max = 90
    eta_step = 0.1
    
    V = np.linspace(0, V_max, int(V_max/V_step))
    eta = np.linspace(0, eta_max, int(eta_max/eta_step))
    x_V, y_eta = np.meshgrid(V, eta)  
    
    Envelope = func_Envelope(mass=mass, eta=y_eta, V=x_V)
    plt.contourf(V, eta, Envelope[0])
    plt.axis('scaled')
    plt.xlabel('Velocity [m/s]')
    plt.ylabel('Rotor Tilt [deg]')
    plt.grid(True)
    plt.show()