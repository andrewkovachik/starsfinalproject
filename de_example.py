import numpy as np
import matplotlib.pyplot as plt
import desolver as de

# Note: This example does not actual need t as an input
#       but is included for consistancy

# Create a new DE solver for a damped oscillating spring
damped_spring = de.DifferentialEquation("Damped Spring")
# Damped spring solution of d^2/dt^2(x) = -b/m d\dt(x) - kx/m
# The lambda is a variable that can be used like a function.
# This defines the variable to have a function which it can
# use to solve the DE.
# It expects that dd is derivatives (0 is x values, 1 is first deriv)
damped_spring.set_derivative_relation(
        lambda dd, t, const: (-const['k']*dd[0]-const['b']*dd[1])/const['m']
        )
# Define characteristic constants
constants = {
        'k': 0.5,
        'b': 0.1,
        'm': 1
        }
# Give initial conditions of outstretched non moving spring
damped_spring.set_boundaries([5, 0])

# Use inital conditions and DE to solve the unknown variable
damped_spring.highest_order_boundary(0, constants)

# Cool we set up the equation lets look at it
print(damped_spring)
#Outputs
# Derivative 0: 5.00	Derivative 1: 0.00	Derivative 2: -2.50

# Sweet now lets take the first step forward!
damped_spring.solve_differential_step(0, 0.1, constants)

# Lets see how it changed!
print(damped_spring)
# Outputs
# Derivative 0: 5.00	Derivative 1: -0.25	Derivative 2: -2.48

# Lets manually step one more time
damped_spring.solve_differential_step(0, 0.1, constants)
print(damped_spring)
# Outputs
# Derivative 0: 4.97	Derivative 1: -0.50	Derivative 2: -2.44

# Everything is changing in a very expected way. Lets step forward in
# constant steps now.

for n in range(48):
    damped_spring.solve_differential_step(0, 0.1, constants)
    print(damped_spring)

# Lastly lets make sure this works for different dampening values

def plot_diff_dampenings(spring, dampening, damp_label):
    time = np.arange(0, 2000)
# No dampening (b^2<4mk)
    constants = {
            'k': 0.5,
            'b': dampening,
            'm': 1
            }
# Give initial conditions of outstretched non moving spring
    spring.set_boundaries([5, 0])
    spring.highest_order_boundary(0, constants)
    for t in time:
        spring.solve_differential_step(0, 0.01, constants)

    plt.plot(time, spring.out_val[0,:-1], label=damp_label)


time = np.arange(0, 2000)
# No dampening (b^2<4mk)
plot_diff_dampenings(damped_spring, 0.1, "Low Dampening")
# critical dampening (b^2=4mk)
plot_diff_dampenings(damped_spring, np.sqrt(2), "Critical Dampening")
# Over dampening (b^2>4mk)
plot_diff_dampenings(damped_spring, 2, "Over Dampening")
plt.legend()

plt.savefig("Dampedspring_example.png")
