"""
Class defining a star and it's various differential equations.
"""
import numpy as np
import desolver as de
import math

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
        info = "STAR: %s\n" % (self.name).upper()
        info += "\n".join("Function {:21}: {}".format(key, str(value))
                          for key, value in self.properties.items()
                          if key not in ["energygen", "radius", "pressure"])

        return info

    def __init__(
            self,
            X=0.55,
            Y=0.43,
            Z=0.02,
            cent_density=162200,
            cent_opticaldepth=0,
            cent_temperature=1.5 * 10**7,
            cent_radii=0.01,  #m
            step_size=0.1,
            error_thresh=1e-5,
            name="Generic Star"):
        """
        Initializes star by deffining the equations that make up
        it's stellar structures, and their differential equations"""

        self.name = name
        self.step_size = step_size

        self.cent_radii = cent_radii
        self.cent_density = cent_density
        self.cent_opticaldepth = cent_opticaldepth
        self.cent_temperature = cent_temperature

        self.X = X
        self.Y = Y
        self.Z = Z
        self.mu = (2 * X + 0.75 * Y + 0.5 * Z)**-1
        self.properties = {
            "opacity": 1,
            "gamma": 5 / 3,
            "pressure": np.array([]),
            "energygen": np.array([]),
            "density": de.RungeKutta("Density"),
            "temperature": de.RungeKutta("Temperature"),
            "mass": de.RungeKutta("Mass"),
            "luminosity": de.RungeKutta("Luminosity"),
            "opticaldepth": de.RungeKutta("Optical Depth"),
            "radius": np.array([cent_radii]),
            "radii": cent_radii
        }

        self.property_list = [
            'opticaldepth', 'luminosity', 'mass', 'temperature', 'density'
        ]
        self.error = [0, 0, 0, 0, 0, 0]
        self.error_thresh = error_thresh

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
            lambda dd, r, state: (4 * np.pi * r**2 * state['density'].now(0) * state["energygen"][-1])
        )

        self.properties['mass'].set_derivative_relation(
            lambda dd, r, state: 4 * np.pi * r**2 * state['density'].now(0))

        self.properties['opticaldepth'].set_derivative_relation(
            lambda dd, r, state: state['opacity'] * state['density'].now(0))

        self.properties['temperature'].set_derivative_relation(
            lambda dd, r, state: -min(3 * state['opacity'] * state['density'].now(0) * state['luminosity'].now(0) / (16 * np.pi * a * C * dd[0]**3 * r**2), (1 - 1 / state['gamma']) * dd[0] * G * state['mass'].now(0) * state['density'].now(0) / (state['pressure'][-1] * r**2))
        )

        self.properties['density'].set_derivative_relation(
            lambda dd, r, state: -(G * self.properties['mass'].now(0) * dd[0] / r**2 + self.properties['pressure_temp_grad'] * self.properties['temperature'].now(1)) / self.properties['pressure_density_grad']
        )

    def setup_boundary_conditions(self):
        """
        Assigns boundary conditions based on their current value.
        Recall that mass, and  luminosity have fixed boundary conditions,
        central pressure is the chosen value, and temperature and optical
        depth are adjusted so to satisfy the surface boundary conditions"""

        cent_energy_pp = 1.07 * 10**-7 * (
            self.cent_density / 10**5) * self.X**2 * (
                self.cent_temperature / 10**6)**4
        cent_energy_cno = 8.24 * 10**-26 * (
            self.cent_density / 10**5) * 0.03 * self.X**2 * (
                self.cent_temperature / 10**6)**19.99

        cent_mass = 4 * np.pi * self.cent_radii**3 * self.cent_density / 3
        cent_lum = 4 * np.pi * self.cent_radii * self.cent_density * (
            cent_energy_pp + cent_energy_cno) / 3
        self.properties['luminosity'].set_boundaries([cent_lum])
        self.properties['mass'].set_boundaries([cent_mass])
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
            (self.properties['density'].now(0) / Mp)**(5 / 3) /
            (5 * Me) + self.properties['density'].now(0) * Kb *
            self.properties['temperature'].now(0) /
            (self.mu * Mp) + a * self.properties['temperature'].now(0)**4 / 3)

        self.properties['pressure_temp_grad'] = (
            self.properties['density'].now(0) * Kb / (self.mu * Mp) +
            4 * a * self.properties['temperature'].now(0)**3 / 3)

        self.properties['pressure_density_grad'] = (
            ((3 * np.pi**2)**(2 / 3) * HBAR**2 *
             (self.properties['density'].now(0) / Mp)**(2 / 3) /
             (3 * Me * Mp)) + Kb * self.properties['temperature'].now(0) /
            (self.mu * Mp))

        k_es = 0.02 * (1 + self.X)
        k_ff = 1 * 10**24 * (self.Z + 0.0001) * (
            self.properties['density'].now(0) / 10**3)**0.7 * (
                self.properties['temperature'].now(0))**-3.5
        k_h = 2.5 * 10**-32 * (self.Z / 0.02) * (
            self.properties['density'].now(0) / 10**3)**0.5 * (
                self.properties['temperature'].now(0))**9

        self.properties['opacity'] = (1 / k_h + 1 / max(k_es, k_ff))**-1

        energy_pp = 1.07 * 10**-7 * (
            self.properties['density'].now(0) / 10**5) * self.X**2 * (
                self.properties['temperature'].now(0) / 10**6)**4
        energy_cno = 8.24 * 10**-26 * (
            self.properties['density'].now(0) / 10**5) * 0.03 * self.X**2 * (
                self.properties['temperature'].now(0) / 10**6)**19.99

        self.properties['energygen'] = np.append(self.properties['energygen'],
                                                 energy_pp + energy_cno)

    def step_de(self):
        """
        Solves the current steps for all the Differential Equations,
        and then steps them forward intime. It is important that the
        step is done after they are all calculated for the cases near
        the central radius. At these values some of the DE's change
        very quickly and so their coupled relationship requires that
        they all use the same input values.
        """
        radius = self.properties['radius'][-1]

        for kutta_const in range(6):

            for item in self.property_list:
                self.properties[item].solve_runge_kutta_const(
                    radius, self.step_size, self.properties, kutta_const)

            for item in self.property_list:
                self.properties[item].use_intermediate()
                self.step_non_de()

        for item in self.property_list:
            self.properties[item].solve_rk_step()
            self.properties[item].solve_de_value(radius, self.step_size,
                                                 self.properties)

        for index, item in enumerate(self.property_list):
            self.error[index] = self.properties[item].error

        if max(self.error) <= self.error_thresh:
            self.properties['radius'] = np.append(
                self.properties['radius'],
                self.properties['radius'][-1] + self.step_size)
            for item in self.property_list:
                self.properties[item].add_differential_step()
                self.properties['radii'] = self.properties['radius'][-1]

            if max(self.error) < 0.3 * self.error_thresh:
                self.adjust_step_size()

        else:
            self.adjust_step_size()
            self.step_de()

    def adjust_step_size(self):
        """
        Uses a relatively quick, and smart way of adjusting the step size.
        Can increase or decrease depending on the threshold of the
        ratio
        """
        self.step_size = self.step_size * (
            self.error_thresh / max(self.error))**.2
        if math.isnan(self.step_size):
            raise ValueError()
