import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt

sigma = 5.6703 * 10**-8


def Lum_error(star):

    radius_surface = star.properties["radius"][-1]
    temperature_surface = star.properties["temperature"].data(0)[-1]
    L_star = star.properties["luminosity"].data(0)[-1]
    return (L_star - 4.0 * np.pi * sigma * (radius_surface**2) *
                     (temperature_surface**4)) / (
                         np.sqrt(4.0 * np.pi * sigma * (radius_surface**2) *
                                 (temperature_surface**4) * L_star))


#####inside of make_star after star.solve, array2D[0] = star.properties['radius']######
"""
Takes in intial guesses for central temperature and pressure
as well as the core type and uses them to create a star and save
a text file.
"""
import stellar_properties as starprop
import Use_Data as data


def make_star(central_temperature, central_density, core_type, name):

    rho_c = central_density
    rho_c_low =  rho_c - 0.1 * rho_c
    rho_c_high = rho_c + 0.1 * rho_c
    step = rho_c * 0.2
    tolerance = 0.001
    i = 1

    while ((rho_c_high - rho_c_low) / 2.0) > tolerance:

        print(central_temperature, rho_c, core_type)
        star = starprop.Star(
            cent_density=float(rho_c),
            cent_temperature=float(central_temperature),
            core=core_type,
            name=name)

        star_low = starprop.Star(
            cent_density=float(rho_c_low),
            cent_temperature=float(central_temperature),
            core=core_type,
            name=name)

        star_high = starprop.Star(
            cent_density=float(rho_c_high),
            cent_temperature=float(central_temperature),
            core=core_type,
            name=name)

        good_solve = star.solve()
        good_solve1 = star_low.solve()
        good_solve2 = star_high.solve()

        print("Low Error: ", (Lum_error(star_low)))
        print("Regular Error: ", (Lum_error(star)))
        print("High Error: ", (Lum_error(star_high)))

        if good_solve == True:
            pass

        else:
            print("Something bad happened")

        if good_solve1 == True:
            pass

        else:
            print("Something bad happened on lower limit")

        if good_solve2 == True:
            pass

        else:
            print("Something bad happened on upper limit")

        if Lum_error(star_low) * Lum_error(star) < 0:
            rho_c_high = rho_c
        elif Lum_error(star_high) * Lum_error(star) < 0:
            rho_c_low = rho_c

        elif abs(Lum_error(star_low)) or abs(Lum_error(star_high)) < 1.0:
            if abs(Lum_error(star_low)) < abs(Lum_error(star_high)):
                rho_c_high = (rho_c + rho_c_low) / 2
                rho_c_low = rho_c_low - (rho_c + rho_c_low) / 2
            else:
                rho_c_low = (rho_c + rho_c_high) / 2
                rho_c_high = rho_c_high - (rho_c + rho_c_high) / 2
            
        else:
            if abs(Lum_error(star_low)) < abs(Lum_error(star_high)):
                rho_c_high = rho_c_low
                rho_c_low = rho_c_high - (0.2 * rho_c)
            if abs(Lum_error(star_high)) < abs(Lum_error(star_low)):
                rho_c_low = rho_c_high
                rho_c_high = rho_c_low + (0.2 * rho_c)
        if i > 30:
            break

        rho_c = (rho_c_low + rho_c_high) / 2.0
        i += 1

    star = starprop.Star(
        cent_density=float(rho_c),
        cent_temperature=float(central_temperature),
        core=core_type,
        name=name)

    star.solve()

    save_variable = [
        'opticaldepth', 'temperature', 'density', 'luminosity', 'mass',
        'opticaldepth_deriv', 'temperature_deriv', 'density_deriv',
        'luminosity_deriv', 'mass_deriv', "k_es", "k_ff", "k_h", "opacity",
        "pressure", "pressure_temp_grad", "pressure_density_grad", "energy_pp",
        "energy_cno", "energy_He", "energy_C", "energygen"
    ]

    print("Saving star:", name)
    deriv = 0
    array2D = [[] for i in range(len(save_variable) + 1)]
    array2D[0] = star.properties['radius']

    for index, variable in enumerate(save_variable):
        if "_deriv" in variable:
            derive = 1
            variable = variable.replace("_deriv", "")
        array2D[index + 1] = star.properties[variable].data(deriv)

    print("Writing star:", name)
    data.array2D2txt(array2D, ["radius"] + save_variable, name)
