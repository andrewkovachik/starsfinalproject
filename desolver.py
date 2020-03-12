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

def temperature_gradient(optical_depth, density=10**-6,
                         opacity=3, half_temp=10**4):
    """
    Output the optical temperature gradient based on the optical depth.

    optical_depth(float, array): data to be used to calculate the temperature grad
    density(float): photosphere density used in calculation kg/m^3
    opacity(float): mean opacity m^2/kg
    half_temp(float): T_e in units of kelvin
    """

    temperature_fourth = 3/4 * half_temp**4*(density*opacity*optical_depth+2/3)
    temperature_grad = np.power(temperature_fourth,1/4)

    return temperature_grad

def excited_fraction(optical_depth, density=10**-6, opacity=3, half_temp=10**4):
    """
    Currently the constants are wrong and this does not operate correctly.

    Takes in a scalling optical depth variable, as well as density, opacity,
    and half temperature. Outputs the fraction of ions at the optical depth

    """

    consts = 4.03
    expons = -1.2*10**5

    temperature = half_temp * np.power(3/4*(opacity*density*optical_depth + 2/3),.25)
    exp_term = np.exp(expons/temperature)
    temp_term = np.power(temperature, 3/2)

    B_term = consts*temp_term*exp_term

    fraction = 4*np.exp(expons/temperature)* ( 1- (B_term)/2*(np.sqrt(1+4/(B_term))-1))

    return fraction

def int_optical_depth(fraction_n2, optical_range, density=10**-6,
                    var_opacity=3.5*10**5, opacity=3, half_temp=10**4):

    step_size = optical_range[1] - optical_range[0]
    integrated_depth = np.cumsum(fraction_n2)*step_size*var_opacity*density + optical_range*density*opacity
    regular_depth = optical_range*density*opacity

    return integrated_depth, regular_depth
    


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

def specific_power(values, temperature=5500, unit='0', frequency=False):
    """
    Calculations the specific power based on blackbody's emission
    per unit area, ster, ln(wavelengt) using Plancks' function 

    wavelength(float): wavelength of units described by units
    temperature(float): temperature of blackbody in kelvin
    units(str): standard SI units that wavelength is described in

    """
    values = si_prefix_multiple(values, unit)

    if frequency:
        planck_func = (
            (2*PLANCK_CONS*values**4)/SPEED_LIGHT**2 
        /(np.exp(PLANCK_CONS*values/(BOLTZ_CONS*temperature)) -1)
        )
    else:
        planck_func = (
            (2*PLANCK_CONS*SPEED_LIGHT**2/values**4) 
        /(np.exp(PLANCK_CONS*SPEED_LIGHT/(values*BOLTZ_CONS*temperature)) -1)
        )

    return planck_func

def black_body_emission(values, temperature=5500, unit='0', frequency=False):
    """
    Calculations the intensity of blackbody's emission per unit wavelength
    using Plancks' function 

    wavelength(float): wavelength of units described by units
    temperature(float): temperature of blackbody in kelvin
    units(str): standard SI units that wavelength is described in

    """
    values = si_prefix_multiple(values, unit)

    if frequency:
        planck_func = (
            (2*PLANCK_CONS*values**3)/SPEED_LIGHT**2 
        /(np.exp(PLANCK_CONS*values/(BOLTZ_CONS*temperature)) -1)
        )
    else:
        planck_func = (
            (2*PLANCK_CONS*SPEED_LIGHT**2/values**5) 
        /(np.exp(PLANCK_CONS*SPEED_LIGHT/(values*BOLTZ_CONS*temperature)) -1)
        )

    return planck_func

def numerical_integration(x_val, y_val):

    dx = x_val[1:] - x_val[0:-1]
    dydx = y_val[:-1]*dx
    integral = np.sum(dydx)
    return integral

class DifferentialEquation:
    """Alows DE's and boundary conditions to be entered
    and can solve the DE numerically for the inputed
    x values"""

    
    def __init__(self, boundary_cond, x_val):
        """
        Sets initial values
        boundary_cond (list): 1 fewer element than the DE order
        x_val (nd.array): x values used to calculate DE
        """
        self.boundaries = boundary_cond + [0]
        self.x_val = x_val

    def set_derivative_relation(self, differential_equation):
        """
        Sets the equation used to solve DE. The lambda must take
        in a list of differential forms first and ten the x values
        
        differential_equation (lambda): should be of the form lambda dd, x: ...
        """
        self.de_relation = differential_equation

        self.highest_order_boundary()

    def highest_order_boundary(self):
        """
        Calculates the highest order boundary for the case that it is not known
        """

        self.boundaries[-1] = self.de_relation(self.boundaries, self.x_val[0])

        self.out_val = np.array(self.boundaries).reshape(-1, 1)

    def solve_differential_step(self, set_of_differentials, x_val):
        """Calculates the next highest order derivative from a set
        of inputed lower order derivatives using the differential
        equation

        set_of_differentials (list): list of derivatves with higher order
                                    derivatives last
        x_val (float): value of x to use at this step
        """
        return self.de_relation(set_of_differentials, x_val)
       
    def solve_differential(self):
        """
        Iterates over the requested x values and finds the output values
        by calculating dx's and adding them to the next step. Uses the 
        DE to solve the highest order differential for the next step"""

        for x_val, x_prev in zip(self.x_val[1:], self.x_val[0:-1]):
                dx = x_val-x_prev

                step = np.array(list(reversed(
                    [0]+[dx*derivative
                        for derivative in reversed(self.out_val[1:,-1])]
                    )))
                step = step+self.out_val[:,-1]
                step[-1] = self.solve_differential_step(step, x_val)
                if step[0] <=0:
                    # This is not general and should be removed in the future
                    # Stops our case from having negative pressure
                    step = np.array([0, 0, 0])

                self.out_val  = np.append(self.out_val, step.reshape(-1,1), axis=1)
