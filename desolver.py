import matplotlib.pyplot as plt
import numpy as np

# Enabling using LaTex for making labels pretty, not really necessary
from matplotlib import rc
rc('text', usetex=True)

PLANCK_CONS = 6.626*10**-34
SPEED_LIGHT = 2.998*10**8
BOLTZ_CONS = 1.381*10**-23

def si_prefix_multiple(value, prefix):
    si_prefixes = {
            'P':10**15,
            'T':10**12,
            'G':10**9,
            'M':10**6,
            'k':10**3,
            'h':10**2,
            'da':10**1,
            '0':10**0,
            'd':10**-1,
            'c':10**-2,
            'm':10**-3,
            'u':10**-6,
            'n':10**-9,
            'p':10**-12,
            'f':10**-15}

    return value * si_prefixes[prefix]

def newtons_method(func, x, dx=0.0001, dylimit=0.000001, iter_limit=100,
                   show_hist=False):
    """
    Newtons method for solving a single variable equation that is equal to zero.

   func(lambda): A single variable lambda function that needs `x` as an input
    x(float): Initial guess to the solution of the function
    dx(float): The small step used to calculate the derivatieve of `func`
    dylimit(float): How close to zero is close enough to end loop
    iter_limit(int): How many iterations before forcing stop
    show_hist(Bool): Print history of steps to screen upon completion

    """
    history = [x]
    for iter in range(0, iter_limit):
        slope = (func(x+dx)-func(x))/dx
        x = x - func(x)/slope
        y_value = func(x)
        history += x

        if np.abs(y_value) <= dylimit:
            if history:
                print(history)
            return x

        elif iter > iter_limit:
            print("Exceded iteration limit")
            return x


def numerical_integration(x_val, y_val):

    dx = x_val[1:] - x_val[0:-1]
    dydx = y_val[:-1]*dx
    integral = np.sum(dydx)
    return integral

class DifferentialEquation:
    """Alows DE's and boundary conditions to be entered
    and can solve the DE numerically for the inputed
    x values"""

    def __str__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("Derivative {}: {:0.4f}".format(n,i) for n, i in enumerate(self.val[:,-1])) 

        return info
    def __repr__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("Derivative {}: {:0.3f}".format(n,i) for n, i in enumerate(self.val[:,-1])) 

        return info
    
    def __init__(self, name="DE Solver"):
        """
        Sets initial values
        boundary_cond (list): 1 fewer element than the DE order
        x_val (nd.array): x values used to calculate DE
        """
        self.name = name

    def set_boundaries(self, boundary_cond):
        """
        Sets the boundary conditions to the given inputs. Boundary 
        conditions are assumed to be for all orders excluding
        the highest order term.
        Args:
            boundary_cond (nd.array): input boundary conditions
        """
        self.boundaries = boundary_cond + [0]
        self.val = np.array(self.boundaries).reshape(-1, 1)
  
    def set_derivative_relation(self, differential_equation):
        """
        Sets the equation used to solve DE. The lambda must take
        in a list of differential forms first and ten the x values
        
        differential_equation (lambda): should be of the form lambda dd, x: ...
        """
        self.de_relation = differential_equation

    def highest_order_boundary(self, x_val, state_vars):
        """
        Calculates the highest order boundary for the case that it is not known
        """

        self.boundaries[-1] = self.de_relation(
                self.boundaries, x_val, state_vars)

        self.val = np.array(self.boundaries).reshape(-1, 1)

    def solve_differential_step(self, x_val, step_size, state_vars):
        """Calculates the next highest order derivative from a set
        of inputed lower order derivatives using the differential
        equation

        x_val (float): value of x to use at this step
        """

        step = np.array(list(reversed(
            [0]+[step_size*derivative
                for derivative in reversed(self.val[1:,-1])]
            )))
        step = step+self.val[:,-1]
        step[-1] = self.de_relation(step, x_val, state_vars)


        self.val  = np.append(self.val, step.reshape(-1,1), axis=1)
