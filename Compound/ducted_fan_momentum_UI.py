import ducted_fan_calc
import matplotlib.pyplot as plt 
import numpy as np

fan_1 = ducted_fan_calc.Ducted_Fan_1(mtow=float(718/4))

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


x = np.arange(0,60, 5)
y = []
for i in x :
    #assuming Cd is constant
    fan_2 = ducted_fan_calc.Ducted_Fan_2(mtow=float(718/4), Cd0=0.25, V= i, related_fan=fan_1)
    # append my y list 
    #temp1 =  fan_2.calc_v_f()
    #temp2 = fan_2.calc_rotor_alpha()
    #temp3 = fan_2.calc_V_horizontal()
    y.append(fan_2.calc_p_idf())
#plot graph 
y = np.array(y)/int(1000)

# 1. Prepare your data
# 2. Create the plot
plt.plot(x, y)  # You can use other types like plt.scatter, plt.bar, etc.

# 3. Customize the plot
# plt.title('My Graph')              # Title of the graph
plt.xlabel('Velocity [m/s]')         # Label for the x-axis
plt.ylabel('Power Required [kW]')         # Label for the y-axis
plt.legend()                       # Add a legend (if needed)
plt.grid(True)                     # Add a grid (optional)

# 4. Display the plot
plt.show()



# fan_3 = ducted_fan_calc.Ducted_Fan_3(mtow=float(718/4), gamma=2, related_fan1=fan_1, related_fan2 = fan_2 )
# print(f"V_c_slow: {fan_3.calc_V_c_slow()}")
# print(f"V_c_fast: {fan_3.calc_V_c_fast()}")
