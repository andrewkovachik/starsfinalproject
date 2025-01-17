import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt

sigma = 5.6703 * 10**-8


def Lum_error(star):

    radius_surface = star.properties["radius"][-1]
    temperature_surface = star.properties["temperature"].data(0)[-1]
    L_star = star.properties["luminosity"].data(0)[-1]
    L_bolt = 4.0 * np.pi * sigma * (radius_surface**2) * (temperature_surface**4)
    ratio = (L_star - L_bolt) / ( np.sqrt(L_bolt * L_star))
    return ratio


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
    rho_c_low =  300
    if core_type == "Hydrogen":
        rho_c_high = 500000
    if core_type == "Helium":
        rho_c_high = 7000000000
    if core_type == "Carbon":
        rho_c_high = 90000000000
    tolerance = 0.0001
    rho_tolerance = 0.000001
    i = 1
    error =10000

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

    low_err = Lum_error(star_low)
    reg_err = Lum_error(star)
    high_err = Lum_error(star_high)

    while abs(error) > tolerance:

        print("Low: ", rho_c_low, low_err)
        print("Med: ", rho_c, reg_err)
        print("Hig: ", rho_c_high, high_err)
        error = reg_err

        if np.abs(reg_err) < tolerance or abs(rho_c_high-rho_c_low) < rho_tolerance:
            print(reg_err<tolerance, abs(rho_c_high-rho_c_low) < rho_tolerance)
            break

        if np.sign(reg_err) == np.sign(low_err):
            rho_c_low = rho_c
            low_err = reg_err

        else:
            rho_c_high = rho_c
            high_err = reg_err

        if i > 60:
            print("Outside of tolerance")
            break

        rho_c = (rho_c_high+rho_c_low)/2
        star = starprop.Star(
           cent_density=float(rho_c),
            cent_temperature=float(central_temperature),
            core=core_type,
            name=name)
        good_solve = star.solve()
        reg_err = Lum_error(star)

        i += 1


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
        else:
            deriv = 0
        array2D[index + 1] = star.properties[variable].data(deriv)

    print("Writing star:", name)
    data.array2D2txt(array2D, ["radius"] + save_variable, name)

    return rho_c
