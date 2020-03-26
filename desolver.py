import matplotlib.pyplot as plt
import numpy as np


def si_prefix_multiple(value, prefix):
    """
    Conversion function for different SI prefixes"""
    si_prefixes = {
        'P': 10**15,
        'T': 10**12,
        'G': 10**9,
        'M': 10**6,
        'k': 10**3,
        'h': 10**2,
        'da': 10**1,
        '0': 10**0,
        'd': 10**-1,
        'c': 10**-2,
        'm': 10**-3,
        'u': 10**-6,
        'n': 10**-9,
        'p': 10**-12,
        'f': 10**-15
    }

    return value * si_prefixes[prefix]


class DifferentialEquation:
    """Alows DE's and boundary conditions to be entered
    and can solve the DE numerically for the inputed
    x values"""

    def __str__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("Derivative {}: {:.2e}".format(n, i)
                         for n, i in enumerate(self.val[:, -1]))

        return info

    def __repr__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("Derivative {}: {:.2e}".format(n, i)
                         for n, i in enumerate(self.val[:, -1]))

        return info

    def __init__(self, name="DE Solver"):
        """
        Sets initial values
        boundary_cond (list): 1 fewer element than the DE order
        x_val (nd.array): x values used to calculate DE
        """
        self.name = name
        self.boundaries = []
        self.val = []
        self.de_relation = None
        self.step = []

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

        Args:
            differential_equation (lambda dd, x, state_var):
        """
        self.de_relation = differential_equation

    def highest_order_boundary(self, x_val, state_vars):
        """
        Calculates the highest order boundary for the case that it is not known

        Args:
            x_val (float): Current x_value or equivilant
            state_vars (list or dict): Callable by a lambda function
        """

        self.boundaries[-1] = self.de_relation(self.boundaries, x_val,
                                               state_vars)

        self.val = np.array(self.boundaries).reshape(-1, 1)

    def runge_kutta(self, step, x_val, step_size, state_vars):
        """
        Runge-kutta method provides a correction factor to a first order
        PDE. A fourth order solution is used here. 

        Args:
            step (nd.array): The previously determined increment of steps
            x_val (float): dependent variable at location of evaluation
            step_size (float): Small step forward being used for calculation
            state_vars (dict): Set of constants that can be used by de_relation
        """
        k1 = step[0]
        k2 = step_size * self.de_relation(self.val[:, -1] + (k1 / 2),
                                          x_val + step_size / 2, state_vars)
        k3 = step_size * self.de_relation(self.val[:, -1] + (k2 / 2),
                                          x_val + step_size / 2, state_vars)
        k4 = step_size * self.de_relation(self.val[:, -1] +
                                          (k3), x_val + step_size, state_vars)

        adjusted_step = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        step[0] = adjusted_step
        return step

    def solve_differential_step(self,
                                x_val,
                                step_size,
                                state_vars,
                                auto_add=True):
        """Calculates the next highest order derivative from a set
        of inputed lower order derivatives using the differential
        equation

        x_val (float): value of x to use at this step
        step_size (float): delta x used to step forward
        state_vars (list or dict): callable by state vars
        auto_add (bool): Whether to immediately add the small step forward,
            it's okay to add immediately if it is an uncoupled DE. If they are
            coupled it should be stepped forward afterwards by itself.
        """

        self.step = np.array(
            list(
                reversed([0] + [
                    step_size * derivative
                    for derivative in reversed(self.val[1:, -1])
                ])))
        # Runge-kutta step has to go here
        if len(self.step) == 2:
            self.step = self.runge_kutta(self.step, x_val, step_size,
                                         state_vars)

        self.step = self.step + self.val[:, -1]
        self.step[-1] = self.de_relation(self.step, x_val, state_vars)

        if auto_add:
            self.add_differential_step()

    def add_differential_step(self):
        """
        Adds the small step as a new set of value to the outvalues
        """
        self.val = np.append(self.val, self.step.reshape(-1, 1), axis=1)
