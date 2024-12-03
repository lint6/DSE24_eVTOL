import numpy as np
import matplotlib.pyplot as plt

## Alexander made this

class PowerAnalysis:        ## IT MUST BE NOTED THAT TAKE OFF AND LANDING ARE NOT INCLUDED -> find time and power required for speed up and speed down the rotors

                            ## ALSO NoTE: all the Velocities are HORIZONTAL VELOCITIES!!!! (this is )
    
    def __init__(self, V_roc = 0.76 , rho = 1.225, eta_p = 0.9, V_climb1 = 30, V_cruise = 40, V_loiter = 6, V_descent = 35, V_climb2 = 32, Volts = 840, motor_eta = 0.85):     #Yunjae, help a homie out brother (Velocities here are just random values to test)

        #values 
        self.Vroc = V_roc       #Vertcial rate of climb
        self.eta_p = eta_p      #propulsion efficiency
        self.rho = rho
        self.Volts = Volts
        self.V_climb1 = V_climb1
        self.V_cruise = V_cruise
        self.V_loiter = V_loiter
        self.V_descent = V_descent
        self.V_climb2 = V_climb2
        self.motor_eta = motor_eta

        #add power required for take off, climb, cruise, descent, loiter, land (corresponding phases can be found in the energy definition)
        self.t1 = None          #[s]
        self.t2 = None
        self.t3 = None
        self.t4 = None
        self.t5 = None
        self.t6 = None
        self.t7 = None
        self.t8 = 1800          # Loiter time, needs to be a set value by the group for max endurance
        self.t9 = None
        self.t10 = None
        self.t11 = None
        self.t12 = None
        self.t13 = None
        self.t14 = None
        self.t_total = None
        self.d_total = None
        self.d_loiter = None    #calculated distance loitered, dependant on V_loiter(endurance) and t8(loiter time)

        ##import power required from somewhere
        self.P1 = 60000          #[W]
        self.P2 = 85000
        self.P3 = 65000
        self.P4 = 100000
        self.P5 = 75000
        self.P6 = 40000
        self.P7 = 65000
        self.P8 = 70000          
        self.P9 = 100000
        self.P10 = 75000
        self.P11 = 40000
        self.P12 = 65000
        self.P13 = 37500
        self.P14 = 60000
        self.P_req_tot = None

        #for the energy
        self.E1 = None          #[s]
        self.E2 = None
        self.E3 = None
        self.E4 = None
        self.E5 = None
        self.E6 = None
        self.E7 = None
        self.E8 = None          # Loiter time, needs to be a set value by the group for max endurance
        self.E9 = None
        self.E10 = None
        self.E11 = None
        self.E12 = None
        self.E13 = None
        self.E14 = None
        self.E_total = None
        


        #NEED TO BE ADDED - import the V_climb1 FROM PERFORMANCE!! - Ask Yunjae!
        #V_climb1 = 36.0       #[m/s] test value


    
    def calculate_missionphase_time(self):
        # HIGE 1
        self.t1 = 15                #[s] - Given
        # Vertical Climb
        delta_h1 = 60.0     #[m] - Given
        self.t2 = delta_h1 / self.Vroc
        # HOGE 1  
        self.t3 = 10             #[s] - Given
        # Steady Climb 1
        h2 = 300        #[m] - Given
        delta_h2 = h2-delta_h1      #[m] height between Hoge1 and cruise 1
        gamma_1 = 9         #[deg] flight path angle
        d_cl1 = delta_h2 / np.tan(gamma_1*np.pi/180)      #Horizontal distance after converting deg to rad
        self.t4 = d_cl1 / self.V_climb1        #time of climb 1 with horizontal climb speed V_climb1
        # Descent 1
        h3 = 30         #[m]       Hoge 2 altitude 
        roc_descent = 7.6      #[m/s]   -    given
        delta_h3 = h2 - h3      #[m] height between cruise 1 and Hoge2
        self.t6 = delta_h3 / roc_descent     #[s]
        d_ds1 = self.V_climb1 * self.t6    # horizontal distance covered during descent
        # Cruise 1
        d_tothalf = 30000       #[m] - Given distance from wright memorial to alligators
        d_cr1 = d_tothalf - d_cl1 - d_ds1       #[m] Horizontal cruise1 distance
        self.t5 = d_cr1 / self.V_cruise       #[s] time needed for cruise1, depending on the cruise speed that has to be implimented!
        # HOGE 2
        self.t7 = 30     #[s]
        # Loitering
        self.t8 = self.t8     #DEPENDS ON MAX ENDURENCE!!-> to be designed for
        self.d_loiter = self.V_loiter * self.t8
        # Steady Climb 2
        delta_h4 = h2 - h3
        gamma_2 = 9     #[deg] flight path angle
        d_cl2 = delta_h4 / np.tan(gamma_2*np.pi/180)        #[m] horizontal distance of climb 2
        self.t9 = d_cl2 / self.V_climb2       #[s] 
        # Descent 2
        delta_h5 = h2-delta_h1
        gamma_3 = 5     #[deg] flight path angle down
        d_ds2 = delta_h5 / np.tan(gamma_3*np.pi/180)        # horizontal distance covered during descent 2
        self.t11 = d_ds2 / self.V_descent     #[s] time to cover horizontal distance according to horizontal descent speed
        # Cruise 2
        d_cr2 = d_tothalf - d_cl2 - d_ds2       #[m] Horizontal cruise2 distance
        self.t10 = d_cr2 / self.V_cruise       #[s] time needed for cruise2, depending on the cruise speed that has to be implimented!
        # HOGE 3
        self.t12 = 10        #[s]
        # Vertical Descent
        Vrod = 0.5      #[m/s]  Vertical rate of descent
        self.t13 = delta_h1 / Vrod       #[s] Vertical descent time
        # HIGE 2
        self.t14 = 15        #[s]

        self.t_total = self.t1 + self.t2 + self.t3 + self.t4 + self.t5 + self.t6 + self.t7 + self.t8 + self.t9 + self.t10 + self.t11 + self.t12 + self.t13 + self.t14
        self.d_total = 2 * d_tothalf + self.d_loiter

        return self.t1, self.t2, self.t3, self.t4, self.t5, self.t6, self.t7, self.t8, self.t9, self.t10, self.t11, self.t12, self.t13, self.t14, self.t_total, self.d_total, self.d_loiter
    
    
    def calculate_energy_required(self):
        self.E1 = self.t1/3600 * self.P1/self.motor_eta                 #Energy required for HIGE 1                 [Wh]
        self.E2 = self.t2/3600 * self.P2/self.motor_eta                 #Energy required for Vertical Climb
        self.E3 = self.t3/3600 * self.P3/self.motor_eta                 #Energy required for HOGE 1
        self.E4 = self.t4/3600 * self.P4/self.motor_eta                 #Energy required for Steady Climb 1
        self.E5 = self.t5/3600 * self.P5/self.motor_eta                 #Energy required for Cruise 1
        self.E6 = self.t6/3600 * self.P6/self.motor_eta                 #Energy required for Descent 1
        self.E7 = self.t7/3600 * self.P7/self.motor_eta                 #Energy required for HOGE 2
        self.E8 = self.t8/3600 * self.P8/self.motor_eta                 #Energy required for Loiter
        self.E9 = self.t9/3600 * self.P9/self.motor_eta                 #Energy required for Steady Climb 1
        self.E10 = self.t10/3600 * self.P10/self.motor_eta                 #Energy required for Cruise 2
        self.E11 = self.t11/3600 * self.P11/self.motor_eta                 #Energy required for Descent 2
        self.E12 = self.t12/3600 * self.P12/self.motor_eta                 #Energy required for HOGE 3
        self.E13 = self.t13/3600 * self.P13/self.motor_eta                 #Energy required for Vertical Descent
        self.E14 = self.t14/3600 * self.P14/self.motor_eta                 #Energy required for HIGE 2
        self.E_total = self.E1 + self.E2 + self.E3 + self.E4 + self.E5 + self.E6 + self.E7 + self.E8 + self.E9 + self.E10 + self.E11 + self.E12 + self.E13 + self.E14

        return self.E1, self.E2, self.E3, self.E4, self.E5, self.E6, self.E7, self.E8, self.E9, self.E10, self.E11, self.E12, self.E13, self.E14, self.E_total
    
    def calculate_amps(self):                           #[A]
        self.A1 = self.P1 / (self.Volts * self.motor_eta)
        self.A2 = self.P2 / (self.Volts * self.motor_eta) 
        self.A3 = self.P3 / (self.Volts * self.motor_eta)
        self.A4 = self.P4 / (self.Volts * self.motor_eta)
        self.A5 = self.P5 / (self.Volts * self.motor_eta)
        self.A6 = self.P6 / (self.Volts * self.motor_eta)
        self.A7 = self.P7 / (self.Volts * self.motor_eta)
        self.A8 = self.P8 / (self.Volts * self.motor_eta)
        self.A9 = self.P9 / (self.Volts * self.motor_eta)
        self.A10 = self.P10 / (self.Volts * self.motor_eta)
        self.A11 = self.P11 / (self.Volts * self.motor_eta)
        self.A12 = self.P12 / (self.Volts * self.motor_eta)
        self.A13 = self.P13 / (self.Volts * self.motor_eta)
        self.A14 = self.P14 / (self.Volts * self.motor_eta)
        self.A_max = max(self.A1, self.A2, self.A3, self.A4, self.A5, self.A6, self.A7, self.A8, self.A9, self.A10, self.A11, self.A12, self.A13, self.A14)
        return self.A1, self.A2, self.A3, self.A4, self.A5, self.A6, self.A7, self.A8, self.A9, self.A10, self.A11, self.A12, self.A13, self.A14, self.A_max
    


    
    # def calculate_total_p_req(self):
        
    #     P_req_tot = self.P1 + self.P2 + self.P3 + self.P4 + self.P5 + self.P6 + self.P7 + self.P8 + self.P9 + self.P10 + self.P11 + self.P12 + self.P13 + self.P14
        

    #     return self.P_req_tot






if __name__ == '__main__':      #printing the P_req vs t plot
    power = PowerAnalysis()
    
    #initiate Power Diagram
    data_P = [power.P1/power.motor_eta, power.P2/power.motor_eta, power.P3/power.motor_eta, power.P4/power.motor_eta, power.P5/power.motor_eta, power.P6/power.motor_eta, power.P7/power.motor_eta, power.P8/power.motor_eta, power.P9/power.motor_eta, power.P10/power.motor_eta, power.P11/power.motor_eta, power.P12/power.motor_eta, power.P13/power.motor_eta, power.P14/power.motor_eta]
    power.calculate_missionphase_time()
    ts = [0, power.t1, power.t2, power.t3, power.t4, power.t5, power.t6, power.t7, power.t8, power.t9, power.t10, power.t11, power.t12, power.t13, power.t14]
    bins = [0]

    for i in range(len(ts)-1):
        bins.append(float(bins[i]) + float(ts[i+1]))
    power_brute = []
    for i in range(len(bins)-1):
        count = 0
        while count < data_P[i]:
            power_brute.append(bins[i])
            count += 1

    fig, axs = plt.subplots(2, 1)  #2 rows, 1 columns of subplots

    #plot first diagram
    axs[0].hist(power_brute, bins=bins, edgecolor='olive', color='darkkhaki')
    axs[0].set_xlabel('Time [s]')
    axs[0].set_ylabel('Power Required [W]')
    axs[0].set_title('Power Required For Each Flight Phase')

    #Height labels Power
    bin_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)]  # Compute bin centers
    for center, count in zip(bin_centers, data_P):
        axs[0].text(center, count + 0.5, f'{int(count)} W', ha='center', va='bottom', fontsize=10, color='red', rotation=90)  # Count above bin

    # Add vertical labels for each bin (Phase Duration)
    bin_centers = [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)]  # Compute bin centers
    for center, time_label in zip(bin_centers, ts[1:]):  # Use time intervals as labels
        axs[0].text(center, -1, f'{round(time_label, 2)} s', ha='center', va='top', fontsize=10, rotation=90)

    # Adjust plot limits to fit vertical labels
    axs[0].set_ylim(bottom=-20000)  # Extend y-axis to make room for vertical labels
    axs[0].set_ylim(top=150000)  # Extend y-axis to make room for vertical labels
    # ------------------------------------------------------------------------------------------------------------------------------------

    #initiate Energy Diagram
    power.calculate_energy_required()
    data_E = [power.E1, power.E2, power.E3, power.E4, power.E5, power.E6, power.E7, power.E8, power.E9, power.E10, power.E11, power.E12, power.E13, power.E14]
    power.calculate_missionphase_time()
    tse = [0, power.t1, power.t2, power.t3, power.t4, power.t5, power.t6, power.t7, power.t8, power.t9, power.t10, power.t11, power.t12, power.t13, power.t14]
    binz = [0]

    for i in range(len(ts)-1):
        binz.append(float(binz[i]) + float(tse[i+1]))
    power_brutez = []
    for i in range(len(binz)-1):
        count = 0
        while count < data_E[i]:
            power_brutez.append(binz[i])
            count += 1
  

    #plot second diagram
    axs[1].hist(power_brutez, bins=binz, edgecolor='black', color='gray')
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Energy Required [Wh]')
    axs[1].set_title('Energy Required For Each Flight Phase')

    #Height labels Energy
    binz_centers = [0.5 * (binz[i] + binz[i + 1]) for i in range(len(binz) - 1)]  # Compute bin centers
    for center, count in zip(binz_centers, data_E):
        axs[1].text(center, count + 0.5, f'{int(count)} Wh', ha='center', va='bottom', fontsize=10, color='red', rotation=90)  # Count above bin

    # Add vertical labels for each bin (Phase Duration)
    bin_centers = [0.5 * (binz[i] + binz[i + 1]) for i in range(len(binz) - 1)]  # Compute bin centers
    for center, time_label in zip(binz_centers, tse[1:]):  # Use time intervals as labels
        axs[1].text(center, -1, f'{round(time_label, 2)} s', ha='center', va='top', fontsize=10, rotation=90)

    # Adjust plot limits to fit vertical labels
    axs[1].set_ylim(bottom=-20000)  # Extend y-axis to make room for vertical labels
    axs[1].set_ylim(top=55000)  # Extend y-axis to make room for vertical labels
    # Show the plot
    plt.show()
     
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #print the energy required per phase
    print('Total Energy Required : ', power.calculate_energy_required()[-1], '[Wh]', 'or', power.calculate_energy_required()[-1]/1000, '[kWh]')

    #print the max amps
    print('Max Ampere: ', power.calculate_amps()[-1], '[A]')
    #print the max amps
    print('Max Voltage: ', power.Volts, '[V]')

