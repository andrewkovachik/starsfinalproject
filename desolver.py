import numpy as np


class DifferentialEquation:
    """Alows DE's and boundary conditions to be entered
    and can solve the DE numerically for the inputed
    x values"""

    def __str__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("{:30} {}: {:.2e}".format(self.name, n, i)
                         for n, i in enumerate(self.val[:, -1]))

        return info

    def __repr__(self):
        """
        Print out useful information for debugging
        """
        info = "\t".join("{:30} {}: {:.2e}".format(self.name, n, i)
                         for n, i in enumerate(self.val[:, -1]))

        return info

    def __init__(self, name="DE Solver"):
        """
        Sets initial values
        """
        self.name = name
        self.boundaries = []
        self.val = []
        self.de_relation = None
        self.step = []
        self.current = []

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
        self.current = self.boundaries

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
        self.current = self.boundaries

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

        self.step = self.step + self.val[:, -1]
        self.step[-1] = self.de_relation(self.step, x_val, state_vars)

        if auto_add:
            self.add_differential_step()

    def add_differential_step(self):
        """
        Adds the small step as a new set of value to the outvalues
        """
        self.val = np.append(self.val, self.step.reshape(-1, 1), axis=1)
        self.current = self.step

    def now(self, order=None):
        """
        Get the value of the Differentiall equation right now
        """

        if order is not None:
            return np.copy(self.current[order])

        else:
            return np.copy(self.current)

    def data(self, order=None):
        """
        Returns full rows of data
        """

        if order:
            return np.copy(self.val[order, :])
        else:
            return np.copy(self.val[0, :])


class RungeKutta(DifferentialEquation):
    def __init__(self, name="DE Solver"):
        """
        Sets initial values
        """
        self.kutta = [0, 0, 0, 0, 0, 0]

        self.y_adj = [
            lambda y, k: y,
            lambda y, k: y + k[0] / 4,
            lambda y, k: y + k[0] * 3 / 32 + k[1] * 9 / 32,
            lambda y, k: y + k[0] * 1932 / 2197 - k[1] * 7200 / 2197 + k[2] * 7296 / 2197,
            lambda y, k: y + k[0] * 439 / 216 - k[1] * 8 + k[2] * 3680 / 513 - k[3] * 845 / 4104,
            lambda y, k: y - k[0] * 8 / 27 + k[1] * 2 - k[2] * 3544 / 2565 + k[3] * 1859 / 4104 - k[4] * 11 / 40,
            lambda y, k: y
        ]
        self.x_adj = [
            lambda x, step: x,
            lambda x, step: x + step / 4,
            lambda x, step: x + step * 3 / 8,
            lambda x, step: x + step * 12 / 13,
            lambda x, step: x + step,
            lambda x, step: x + step / 2
        ]

        self.intermediate = []
        self.hold = []
        self.derivative_hold = 0
        self.error = 0
        super().__init__(name)

    def solve_runge_kutta_const(self, x_val, step_size, state_vars,
                                kutta_const):
        """
        Runge-kutta method provides a correction factor to a first order
        PDE. A fourth order solution is used here. 

        Args:
            x_val (float): dependent variable at location of evaluation
            step_size (float): Small step forward being used for calculation
            state_vars (dict): Set of constants that can be used by de_relation
            kutta_const (int): Kutta constant is being solved for

        Returns:
            (nd.array): Adjusted step after making runge-kutta correction
        """
        if kutta_const == 0:
            self.hold = self.now()
            self.intermediate = np.copy(self.hold)

        x_adj = self.x_adj[kutta_const](x_val, step_size)

        result = self.de_relation(self.intermediate, x_adj, state_vars)
        self.intermediate[1] = result

        self.kutta[kutta_const] = step_size * self.intermediate[1]

        self.intermediate[0] = self.y_adj[kutta_const + 1](self.hold[0],
                                                           self.kutta)

    def use_intermediate(self):
        """
        Sets the current value of the DE to be that which is calculated from 
        Runge-kutta method. Usefull for when solving many DE's with Runge-kutta
        """
        self.current = np.copy(self.intermediate)

    def use_original(self):
        """
        Sets the current value to return back to the value being held in the
        hold variable
        """
        self.current = np.copy(self.hold)

    def solve_rk_step(self):
        """
        Use the calculated runge-kutta constants to determine the 4th and 5th
        order solutions as well as their error between the two.
        """

        self.kutta_5th_sol = (
            self.hold[0] + self.kutta[0] * 16 / 135 +
            self.kutta[2] * 6656 / 12825 + self.kutta[3] * 28561 / 56430 -
            self.kutta[4] * 9 / 50 + self.kutta[5] * 2 / 55)

        self.kutta_4th_sol = (
            self.hold[0] + self.kutta[0] * 25 / 216 +
            self.kutta[2] * 1408 / 2565 + self.kutta[3] * 2197 / 4104 -
            self.kutta[4] / 5 + self.kutta[5] * 0)

        self.step = np.array([self.kutta_4th_sol, 0])

        self.error = abs(
            (self.kutta_4th_sol - self.kutta_5th_sol) / self.kutta_5th_sol)

    def solve_de_value(self, x_val, step_size, state_vars):
        """
        Adjust DE value at the current point after the kutta-de step.
        Needed because some DE's will reference other DE value and so
        we cannot wait for them to self update
        """

        self.step[1] = self.de_relation(self.step, x_val + step_size,
                                        state_vars)
