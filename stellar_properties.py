"""
Class defining a star and it's various differential equations.
"""
import numpy as np
import desolver as de

# m/s^2
C = 2.98 * 10**8
#m^3/kg/s^2
G = 6.6741 * 10**-11
# m^2*kg/s
HBAR = 1.055 * 10**-34

Me = 9.11 * 10**-31
Mp = 1.67 * 10**-27
#Boltzman J/K
Kb = 1.381 * 10**-23
# FIND THIS
# J/k^4/m^3
a = 7.566 * 10**-16


class Star:
    """
    Class definining star. Can calculate many different
    properties of the star at the various lengths across
    it's radius
    """

    def __str__(self):
        """
        Print out useful information for debugging
        """
        info = "STAR: %s\n"%(self.name).upper()
        info += "\n".join("Function {:21}: {}".format(key, str(value))
                         for key, value in self.properties.items() if
                         key not in ["pressure", "radius", "energygen"])

        return info

    def __init__(self,
                 X=0.55,
                 Y=0.43,
                 Z=0.02,
                 cent_density=162200,
                 cent_opticaldepth=2 / 3,
                 cent_temperature=1.5*10**7,
                 step_size=0.1,
                 name="Generic Star"):
        """
        Initializes star by deffining the equations that make up
        it's stellar structures, and their differential equations"""

        self.name = name
        self.step_size = step_size

        self.cent_density = cent_density
        self.cent_opticaldepth = cent_opticaldepth
        self.cent_temperature = cent_temperature

        self.X = X
        self.Y = Y
        self.Z = Z
        self.mu = (2 * X + 0.75 * Y + 0.5 * Z)**-1
        self.properties = {
            "opacity": 1,
            "gamma": 5/3,
            "pressure": np.array([]),
            "energygen": np.array([]),
            "density": de.DifferentialEquation("Density"),
            "temperature": de.DifferentialEquation("Temperature"),
            "mass": de.DifferentialEquation("Mass"),
            "luminosity": de.DifferentialEquation("Luminosity"),
            "opticaldepth": de.DifferentialEquation("Optical Depth"),
            "radius": np.array([0.001])
        }

        self.setup_stellar_equations()
        self.setup_boundary_conditions()
        self.step_non_de()

    def setup_stellar_equations(self):
        """
        Assigns the stellar properties their differential equation.
        Equations are taken from the assignment manual. It is assumed
        that dd is the derivatives of that stellar property, r is radius
        and state includes a dictionary of all other stellar properties. 
        In the case that the property is another DE variable, the first
        index defines which derivative is grabbed, and the second index 
        defines the which element to grab
        """
        self.properties['luminosity'].set_derivative_relation(
            lambda dd, r, state: 4 * np.pi * r**2 * state['density'].val[0, -1]* state["energygen"][-1])

        self.properties['mass'].set_derivative_relation(
            lambda dd, r, state: 4 * np.pi * r**2 * state['density'].val[0, -1]
        )

        self.properties['opticaldepth'].set_derivative_relation(
            lambda dd, r, state: state['opacity'] * state['density'].val[0, -1]
        )

        self.properties['temperature'].set_derivative_relation(
            lambda dd, r, state: -min(
                3 * state['opacity'] * state['density'].val[0, -1] * state['luminosity'].val[0, -1] / (16 * np.pi * a * C * dd[0]**3 * r**2),
                (1 - 1 / state['gamma']) * dd[0] * G * state['mass'].val[0, -1] * state['density'].val[0, -1] / (state['pressure'][-1] * r**2))
        )

        self.properties['density'].set_derivative_relation(
            lambda dd, r, state: 0-(G * state['mass'].val[0, -1] * dd[0] / r**2 + state['pressure_temp_grad'] * state['temperature'].val[1, -1]) / state['pressure_density_grad']
        )

    def setup_boundary_conditions(self):
        """
        Assigns boundary conditions based on their current value.
        Recall that mass, and  luminosity have fixed boundary conditions,
        central pressure is the chosen value, and temperature and optical
        depth are adjusted so to satisfy the surface boundary conditions"""

        self.properties['luminosity'].set_boundaries([0])
        self.properties['mass'].set_boundaries([0])
        self.properties['opticaldepth'].set_boundaries(
            [self.cent_opticaldepth])
        self.properties['temperature'].set_boundaries([self.cent_temperature])
        self.properties['density'].set_boundaries([self.cent_density])

    def step_non_de(self):
        """Updates the variables and arrays of variables that do not depend
        on differential equations and can be calculated dirrectly.
        """

        self.properties['pressure'] = np.append(
            self.properties['pressure'], (3 * np.pi**2)**(2 / 3) * HBAR**2 *
            (self.properties['density'].val[0, -1] / Mp)**(5 / 3) /
            (5 * Me) + self.properties['density'].val[0, -1] * Kb *
            self.properties['temperature'].val[0, -1] / (self.mu * Mp) +
            a * self.properties['temperature'].val[0, -1]**4 / 3)

        self.properties[ 'pressure_temp_grad'] = (
                self.properties['density'].val[0, -1] * Kb / ( self.mu * Mp)
                + 4 * a * self.properties['temperature'].val[0, -1]**3)

        self.properties['pressure_density_grad'] = (
                (3 * np.pi**2)**( 2 / 3) * HBAR**2
                * (self.properties['density'].val[0, -1] / Mp)**(2 / 3)
                / (3 * Me * Mp) + Kb * self.properties['temperature'].val[0, -1] / (
                self.mu * Mp))

        k_es = 0.02 * (1 + self.X)
        k_ff = 1 * 10**24 * (self.Z + 0.0001) * (
            self.properties['density'].val[0, -1] / 10**3)**0.7 * (
                self.properties['temperature'].val[0, -1])**-3.5
        k_H = 2.5 * 10**-32 * (self.Z / 0.02) * (
            self.properties['density'].val[0, -1] / 10**3)**0.5 * (
                self.properties['temperature'].val[0, -1])**9
        
        self.properties['opacity'] = (1 / k_H + 1 / max(k_es, k_ff))**-1

        energy_pp = 1.07 * 10**-6 * (
            self.properties['density'].val[0, -1] / 10**5) * self.X**2 * (
                self.properties['temperature'].val[0, -1] / 10**6)**4
        energy_cno = 8.24 * 10**-26 * (
            self.properties['density'].val[0, -1] / 10**5
        ) * 0.03 * self.X**2 * (
            self.properties['temperature'].val[0, -1] / 10**6)**19.99

        self.properties['energygen'] = np.append(self.properties['energygen'],
                                                 energy_pp + energy_cno)

    def step_de(self):
        self.properties['radius'] = np.append(
            self.properties['radius'],
            self.properties['radius'][-1] + self.step_size)
        radius = self.properties['radius'][-1]

        self.properties['opticaldepth'].solve_differential_step(
            radius, self.step_size, self.properties)

        self.properties['luminosity'].solve_differential_step(
            radius, self.step_size, self.properties)

        self.properties['mass'].solve_differential_step(
            radius, self.step_size, self.properties)

        self.properties['temperature'].solve_differential_step(
            radius, self.step_size, self.properties)

        self.properties['density'].solve_differential_step(
            radius, self.step_size, self.properties)

        self.properties['opticaldepth'].solve_differential_step(
            radius, self.step_size, self.properties)
