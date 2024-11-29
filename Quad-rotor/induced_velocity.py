import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants (Example values, adjust as needed)
v_i_hover = 10  # Hover induced velocity (example value)
V_cruise_values = np.linspace(5, 50, 100)  # Range of V_cruise values to evaluate
k = 1.15  # Induced power factor (example value)
T = 718.89 * 9.81 *  (1 + ((0.04 * 718.89) / 718.89))  # Thrust (in N, example value)

# Equation to solve
def equation(v_i_bar, V_cruise, v_i_hover):
    return v_i_bar**4 + (V_cruise / v_i_hover)**2 * v_i_bar**2 - 1

# Arrays to store solutions
v_i_bar_solutions = []
v_i_fl_solutions = []
P_i_solutions = []

# Solve for each V_cruise
for V_cruise in V_cruise_values:
    # Initial guess for v_i_bar (can be adjusted as needed)
    initial_guess = 1.0
    
    # Solve the equation using fsolve
    v_i_bar_solution = fsolve(equation, initial_guess, args=(V_cruise, v_i_hover))
    
    # Append the first (realistic) solution
    v_i_bar = v_i_bar_solution[0]
    v_i_bar_solutions.append(v_i_bar)
    
    # Compute v_i_fl and store
    v_i_fl = v_i_bar * v_i_hover
    v_i_fl_solutions.append(v_i_fl)
    
    # Compute induced power and store
    P_i = k * T * v_i_fl
    P_i_solutions.append(P_i/1000)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(V_cruise_values, P_i_solutions, label=r'$P_{i}$')
plt.xlabel(r'$V_{cruise}$ (m/s)')
plt.ylabel(r'$P_{i}$ (kW)')
plt.title(r'$\mathbf{Induced\ Power\ (P_{i})}$ vs $\mathbf{V_{cruise}}$')
plt.legend()
plt.grid(True)
plt.show()
