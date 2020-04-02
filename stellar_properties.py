"""
Class defining a star and it's various differential equations.
"""
import numpy as np
import math
import desolver as de
import regular_equation as re

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
sigma = 5.67e-8  # W/m^2 * K^-4


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
            X=0.70,
            Y=0.28,
            Z=0.02,
            Xc=0.004,
            cent_density=162200,
            cent_opticaldepth=0,
            cent_temperature=1.5 * 10**7,
            cent_radii=0.01,  #m
            step_size=0.1,
            error_thresh=1e-5,
            max_step=100000,
            min_step=0.001,
            core="Hydrogen",
            #core is one of "Hydrogen", "Helium", "Carbon"
            name="Generic Star"):
        """
        Initializes star by deffining the equations that make up
        it's stellar structures, and their differential equations"""

        self.name = name
        self.step_size = step_size
        self.max_step = max_step
        self.min_step = min_step

        self.cent_radii = cent_radii
        self.cent_density = cent_density
        self.cent_opticaldepth = cent_opticaldepth
        self.cent_temperature = cent_temperature

        self.Xc = Xc
        self.X = X
        self.Y = Y
        self.Z = Z
        self.mu = (2 * X + 0.75 * Y + 0.5 * Z)**-1
        self.core = core
        self.properties = {
            "opacity": re.Equation("Opacity"),
            "k_es": re.Equation("Electron Scattering Opacity"),
            "k_ff": re.Equation("Free Free Scattering Opacity"),
            "k_h": re.Equation("Hydrogen Opacity"),
            "pressure": re.Equation("Pressure"),
            "pressure_temp_grad": re.Equation("Pressure Temp Gradient"),
            "pressure_density_grad": re.Equation("Pressure Density Gradient"),
            "energygen": re.Equation("Energy Generation"),
            "energy_pp": re.Equation("Proton-Proton Energy Generation"),
            "energy_cno": re.Equation("CNO Cycle Energy Generation"),
            "energy_He": re.Equation("Helium Energy Generation"),
            "energy_C": re.Equation("CNO Cycle Energy Generation"),
            "density": de.RungeKutta("Density"),
            "temperature": de.RungeKutta("Temperature"),
            "mass": de.RungeKutta("Mass"),
            "luminosity": de.RungeKutta("Luminosity"),
            "opticaldepth": de.RungeKutta("Optical Depth"),
            "radius": np.array([cent_radii]),
            "radii": cent_radii,
            "gamma": 5 / 3,
        }

        self.de_list = [
            'opticaldepth', 'temperature', 'density', 'luminosity', 'mass'
        ]
        self.eq_list = [
            "k_es", "k_ff", "k_h", "opacity", "pressure", "pressure_temp_grad",
            "pressure_density_grad", "energy_pp", "energy_cno", "energy_He",
            "energy_C", "energygen"
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
            lambda dd, r, state: (4 * np.pi * r**2 * state['density'].now(0) * state["energygen"].now())
        )

        self.properties['mass'].set_derivative_relation(
            lambda dd, r, state: 4 * np.pi * r**2 * state['density'].now(0))

        self.properties['opticaldepth'].set_derivative_relation(
            lambda dd, r, state: state['opacity'].now() * state['density'].now(0)
        )

        self.properties['temperature'].set_derivative_relation(
            lambda dd, r, state: -min(
                3 * state['opacity'].now() * state['density'].now(0) * state['luminosity'].now(0) / (64*np.pi*r**2*sigma * state['temperature'].now(0)**3),
                (1 - 1 / state['gamma']) * state['temperature'].now(0) * G * state['mass'].now(0) * state['density'].now(0) / (state['pressure'].now() * r**2))
        )

        self.properties['density'].set_derivative_relation(
            lambda dd, r, state: -(G * state['mass'].now(0) * state['density'].now(0) / r**2 + state['pressure_temp_grad'].now() * state['temperature'].now(1)) / state['pressure_density_grad'].now()
        )

        self.properties['pressure'].set_equation(
            lambda state: ((3 * np.pi**2)**(2 / 3) * HBAR**2 * (state['density'].now(0) / Mp)**(5 / 3) / (5 * Me) + state['density'].now(0) * Kb * state['temperature'].now(0) / (self.mu * Mp) + a * state['temperature'].now(0)**4 / 3)
        )

        self.properties['pressure_density_grad'].set_equation(
            lambda state: ((3 * np.pi**2)**(2 / 3) * HBAR**2 * (state['density'].now(0) / Mp)**(2 / 3) / (3 * Me * Mp) + Kb * state['temperature'].now(0) / (self.mu * Mp))
        )

        self.properties['pressure_temp_grad'].set_equation(
            lambda state: (state['density'].now(0) * Kb / (self.mu * Mp) + 4 * a * state['temperature'].now(0)**3 / 3)
        )

        self.properties['energy_pp'].set_equation(
            lambda state: (1.07e-7 * (state['density'].now(0) / 1e5) * self.X**2 * (state['temperature'].now(0) / 1e6)**4)
        )

        self.properties['energy_cno'].set_equation(
            lambda state: (8.24e-26 * (state['density'].now(0) / 1e5) * 0.03 * self.X**2 * (state['temperature'].now(0) / 1e6)**19.9)
        )

        self.properties['energy_He'].set_equation(
            lambda state: 3.85e-8 * (state['density'].now(0) / 1e5)**2 * self.Y**3 * (state['temperature'].now(0) / 1e8)**44
        )

        self.properties['energy_C'].set_equation(
            lambda state: 5.0e4 * (state['density'].now(0) / 1e5) * self.Xc**2 * (state['temperature'].now(0) / 1e9)**30
        )

        if self.core == "Hydrogen":
            self.properties['energygen'].set_equation(
                lambda state: state['energy_pp'].now(0) + state['energy_cno'].now(0)
            )

        if self.core == "Helium":
            self.properties['energygen'].set_equation(
                lambda state: state['energy_He'].now(0))

        if self.core == "Carbon":
            self.properties['energygen'].set_equation(
                lambda state: state['energy_C'].now(0))

        self.properties['k_es'].set_equation(lambda state: 0.02 * (1 + self.X))

        self.properties['k_ff'].set_equation(
            lambda state: 1e24 * (self.Z + 0.0001) * (state['density'].now(0) / 1e3)**0.7 * (state['temperature'].now(0))**-3.5
        )

        self.properties['k_h'].set_equation(
            lambda state: 2.5e-32 * (self.Z / 0.02) * (state['density'].now(0) / 1e3)**0.5 * (state['temperature'].now(0))**9
        )

        self.properties['opacity'].set_equation(
            lambda state: (state['k_h'].now()*max(state['k_es'].now(),state['k_ff'].now())/(state['k_h'].now()+max(state['k_es'].now(),state['k_ff'].now()))))

            # (1/state['k_h'].now() + 1/max(state['k_es'].now(), state['k_ff'].now()))**-1 )

    def setup_boundary_conditions(self):
        """
        Assigns boundary conditions based on their current value.
        Recall that mass, and  luminosity have fixed boundary conditions,
        central pressure is the chosen value, and temperature and optical
        depth are adjusted so to satisfy the surface boundary conditions"""
        
        # Ensures that the read values are numbers not strings
        self.cent_density = float(self.cent_density)
        self.cent_temperature = float(self.cent_temperature)

        cent_energy_pp = 1.07 * 10**-7 * (
            self.cent_density / 10**5) * self.X**2 * (
                self.cent_temperature / 10**6)**4
        cent_energy_cno = 8.24 * 10**-26 * (
            self.cent_density / 10**5) * 0.03 * self.X**2 * (
                self.cent_temperature / 10**6)**19.99
        cent_energy_He = 3.85e-8 * (
            self.cent_density / 10**5)**2 * self.Y**3 * (
                self.cent_temperature / 10**8)**44
        cent_energy_C = 5.0e4 * (self.cent_density / 10**5) * self.Xc**2 * (
            self.cent_temperature / 10**9)**30

        cent_mass = 4 * np.pi * self.cent_radii**3 * self.cent_density / 3

        if self.core == "Hydrogen":
            cent_lum = 4 * np.pi * self.cent_radii * self.cent_density * (
                cent_energy_pp + cent_energy_cno) / 3
        if self.core == "Helium":
            cent_lum = 4 * np.pi * self.cent_radii * self.cent_density * (
                cent_energy_He) / 3
        if self.core == "Carbon":
            cent_lum = 4 * np.pi * self.cent_radii * self.cent_density * (
                cent_energy_C) / 3

        self.properties['luminosity'].set_boundaries([cent_lum])
        self.properties['mass'].set_boundaries([cent_mass])
        self.properties['opticaldepth'].set_boundaries(
            [self.cent_opticaldepth])
        self.properties['temperature'].set_boundaries([self.cent_temperature])
        self.properties['density'].set_boundaries([self.cent_density])

    def step_non_de(self, auto_add=True):
        """Updates the variables and arrays of variables that do not depend
        on differential equations and can be calculated dirrectly.
        """

        for equation in self.eq_list:
            self.properties[equation].solve_step(
                self.properties, auto_add=auto_add)

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
        nan_problem = False

        for kutta_const in range(6):

            for item in self.de_list:
                self.properties[item].solve_runge_kutta_const(
                    radius, self.step_size, self.properties, kutta_const)

            self.de_use_intermediate()
            self.step_non_de(auto_add=False)
            self.eq_use_intermediate()

        for item in self.de_list:
            self.properties[item].solve_rk_step()
            self.properties[item].solve_de_value(radius, self.step_size,
                                                 self.properties)

        for index, item in enumerate(self.de_list):
            self.error[index] = self.properties[item].error


        if max(self.error) <= self.error_thresh:

            self.properties['radius'] = np.append(
                self.properties['radius'],
                self.properties['radius'][-1] + self.step_size)
            for item in self.de_list:
                self.properties[item].add_differential_step()
                self.properties['radii'] = self.properties['radius'][-1]

            self.step_non_de(auto_add=True)
            if max(self.error) < 0.1 * self.error_thresh:
                self.adjust_step_size()

        else:
            self.adjust_step_size()
            for item in self.de_list:
                self.properties[item].use_original()
            for item in self.eq_list:
                self.properties[item].use_original()

    def solve(self):
        """
        Runs a loop within itself until it is satisfied with the
        outer-layer
        """

        self.check_stop()

        while self.run:
            self.step_de()
            self.check_stop()
            if len(self.properties['radius']) % 2000 == 0:
                pass
                #print(self)
                #print(self.name, self.dtau)

        self.remove_extra()
        return self.success

    def remove_extra(self):

        tau_infinity = self.properties['opticaldepth'].now(0)

        tau_adjusted = abs((tau_infinity - self.properties['opticaldepth'].data(0)) - (2/3))
        radius_index = np.argmin(tau_adjusted)
        total_length = len(tau_adjusted)

        for item in self.eq_list:
            self.properties[item].val = self.properties[item].val[:radius_index]

        for item in self.de_list:
            self.properties[item].val = self.properties[item].val[:,:radius_index]

        self.properties['radius'] = self.properties['radius'][:radius_index]

    def adjust_step_size(self):
        """
        Uses a relatively quick, and smart way of adjusting the step size.
        Can increase or decrease depending on the threshold of the
        ratio
        """

        if max(self.error) == 0:
            self.step_size = self.step_size*10

        else:
            self.step_size = max(
                self.min_step,
                    min(
                        self.step_size * 0.8 *
                        (self.error_thresh / max(self.error))**.2,
                        self.max_step)
            )


    def de_use_intermediate(self):
        """
        Forces non de equations to return their intermediate step
        """

        for differential in self.de_list:
            self.properties[differential].use_intermediate()

    def eq_use_intermediate(self):
        """
        Forces non de equations to return their intermediate step
        """

        for equation in self.eq_list:
            self.properties[equation].use_intermediate()

    def check_stop(self):
        """
        Checks closeness to tau infinity
        """

        self.dtau = (self.properties['opacity'].now() *
                     (self.properties['density'].now(0))**2 / abs(
                         self.properties['density'].now(1)))

        if self.dtau < 0.00001:
            self.run = False
            self.success = True

        elif len(self.properties['radius']) > 5000:
            print("Stopping based on large number of iterations > 30000")
            self.run = False
            self.success = False

        else:
            self.run = True
