import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy import symbols, diff, lambdify

# Define the function symbolically
variable = symbols('x')
function = variable**2 + variable + 1  # Example: Replace this with your desired function
f = lambdify(variable, function)  # Convert to a numerical function
f_derivative_nd = diff(function, variable)  # Compute the derivative symbolically
f_derivative = lambdify(variable, f_derivative_nd)  # Convert the derivative to a numerical function

# Tangency condition: f'(x) = f(x) / x
def tangent(x):
    if x == 0:
        return np.inf  # Avoid division by zero
    return f_derivative(x) - f(x) / x

# Find the tangency point numerically
initial_guess = 1.0  # Choose a reasonable initial guess
x_tangent = fsolve(tangent, initial_guess)[0]
y_tangent = f(x_tangent)

# Tangent line equation
def tangent_line(x):
    slope = y_tangent / x_tangent
    return slope * x

# Plotting
x_values = np.linspace(0, 10, 500)  # Adjust the range as needed
y_values = f(x_values)

plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values, label='f(x)')
plt.plot(x_values, tangent_line(x_values), '--', label='Tangent Line')
plt.scatter([x_tangent], [y_tangent], color='red', label=f'Tangent Point ({x_tangent:.2f}, {y_tangent:.2f})')

plt.xlabel('x')
plt.ylabel('f(x)')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.title('Function and Tangent Line from Origin')
plt.legend()
plt.grid()
plt.show()

