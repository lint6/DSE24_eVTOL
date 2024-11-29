import ducted_fan_calc

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

fan_2 = ducted_fan_calc.Ducted_Fan_2(mtow=float(718/4), related_fan=fan_1)
print(f"Forward flight induced velocity from fan 2: {fan_2.calc_v_f()}")
print(f"rotor disc angle  of attack from fan 2: {fan_2.calc_rotor_alpha()}")
print(f"Horizontal velocity: {fan_2.calc_V_horizontal()}")



fan_3 = ducted_fan_calc.Ducted_Fan_3(mtow=float(718/4), gamma=3, related_fan1=fan_1, related_fan2 = fan_2 )
print(f"V_c_slow: {fan_3.calc_V_c_slow()}")
print(f"V_c_slow: {fan_3.calc_V_c_fast()}")
