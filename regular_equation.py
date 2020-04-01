import numpy as np


class Equation:
    """
    Allows for entering equations with some useful methods of
    grabbing and changing the values
    """

    def __str__(self):
        """
        Print out useful information for debugging
        """
        info = ("Equation for {:30}: {:.2e}".format(self.name, self.val[-1]))

        return info

    def __repr__(self):
        """
        Print out useful information for debugging
        """
        info = ("Equation for {:30}: {:.2e}".format(self.name, self.val[-1]))

        return info

    def __init__(self, name="Equation"):
        """
        Sets initial values
        """
        self.name = name
        self.val = np.array([])
        self.step = []
        self.current = []

    def set_equation(self, equation):
        """
        Sets the equation used to solve DE. The lambda must take
        in a the x/t value first and then the state variable

        Args:
            differential_equation (lambda state_var: ... )
        """
        self.equation = equation

    def solve_step(self, state, auto_add=True):
        """
        Solves for the current value based on the current x/t value
        and other variables represented by the state variable. Both
        will be handed to the lambda equation in the order of
        lambda self, x_val, state:

        Args:
            state (float, dict, list): Anything callable in the lambda function
            auto_add (bool): Whether to immediately add to the main array
        """

        self.step = self.equation(state)

        if auto_add:
            self.add_step()

        else:
            self.intermediate = np.copy(self.step)

    def add_step(self):
        """
        Adds the small step as a new set of value to the outvalues
        """
        self.val = np.append(self.val, self.step)
        self.current = np.copy(self.step)

    def now(self, order=None):
        """
        Get the value of the equation right now.
        Order does nothing at the moment other than make it match in 
        syntax to desolver
        """

        return np.copy(self.current)

    def use_intermediate(self):
        """
        Sets the current value of the DE to be that which is calculated from 
        Runge-kutta method. Usefull for when solving many DE's with Runge-kutta
        """
        self.hold = np.copy(self.current)
        self.current = np.copy(self.intermediate)

    def use_original(self):
        """
        Sets the current value to return back to the value being held in the
        hold variable
        """
        self.current = np.copy(self.hold)

    def data(self, order=None):
        """
        Returns full rows of data
        """

        if order:
            pass

        return np.copy(self.val[:])
