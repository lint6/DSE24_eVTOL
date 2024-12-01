import ducted_fan_calc
import matplotlib.pyplot as plt 
import numpy as np

def func_min_locator(list1, list2):
    min_list2 = np.min(list2)
    min_list2_index = list(list2).index(min_list2)
    return list1[min_list2_index], min_list2

fan_1 = ducted_fan_calc.Ducted_Fan_1(mass=float(718/4))

# Calculate hover induced velocity
print(f"Hover Induced Velocity: {fan_1.calc_hover_induced_velovity()}")

# Calculate thrust for hover
print(f"Thrust for Hover: {fan_1.calc_thrust_hover()}")

# Calculate power for steady climb
print(f"Power for Steady Climb: {fan_1.calc_power_ideal_steady_climb()}")

print(f"Rate of Climb rate: {fan_1.calc_V_c_kappa()}")

print(f"calc_power_ideal_hover: {fan_1.calc_power_ideal_hover()}")

print(f"calc_disc_loading: {fan_1.calc_disc_loading()}")

print(f"calc_pidd: {fan_1.calc_p_idd()}")


vel = np.arange(0,100,0.01)
power = []
for i in vel :
    #assuming Cd is constant
    fan_2 = ducted_fan_calc.Ducted_Fan_2(mass=float(718/4), Cd0=0.05, V=i, related_fan=fan_1, radius = 0.625)
    power.append(fan_2.calc_p_idf())
#plot graph 
power = np.array(power)/int(1000)

# 1. Prepare your data
# Power lines
grey = (0.2,0.2,0.2)
plt.plot(vel, power.T[0], "-", label="Total")  # You can use other types like plt.scatter, plt.bar, etc.
plt.plot(vel, power.T[1], "--",label="Induced")
plt.plot(vel, power.T[2], "--",label="parasitic")
#Critical Points
vel_power_min, power_min = func_min_locator(vel, power.T[0])
plt.plot(vel_power_min, power_min, '.', label=f'Min. Power \nV = {vel_power_min:.1f}m/s \nP = {power_min:.2f}kW')
plt.plot(vel[0], power[0][0], '.', label=f'Hover Power P = {power[0][0]:.2f}kW')

# 3. Customize the plot
plt.title('Forward Flight Power (Single rotor)')              # Title of the graph
plt.xlabel('Velocity [m/s]')         # Label for the x-axis
plt.ylabel('Power Required [kW]')         # Label for the y-axis
plt.legend()                       # Add a legend (if needed)
plt.grid(True)                     # Add a grid (optional)

# 4. Display the plot
plt.show()



# fan_3 = ducted_fan_calc.Ducted_Fan_3(mtow=float(718/4), gamma=2, related_fan1=fan_1, related_fan2 = fan_2 )
# print(f"V_c_slow: {fan_3.calc_V_c_slow()}")
# print(f"V_c_fast: {fan_3.calc_V_c_fast()}")
