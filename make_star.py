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
    rho_c_low =  rho_c - 0.6 * rho_c
    rho_c_high = rho_c + 0.6 * rho_c
    step = rho_c * 0.2
    tolerance = 0.01
    i = 1
    error =10000

    while error > tolerance:

        print(central_temperature, rho_c_low, core_type)
        print(central_temperature, rho_c, core_type)
        print(central_temperature, rho_c_high, core_type)
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
        all_err = [low_err, reg_err, high_err]

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
            rho_c_low = (rho_c + rho_c_low)/2
        if good_solve2 == True:
            pass
        else:
            print("Something bad happened on upper limit")
            rho_c_high = (rho_c + rho_c_high)/2

        error = Lum_error(star)

        if all(err > 0 for err in all_err):
            diff = rho_c_high - rho_c

            if low_err < reg_err:
                rho_c = rho_c_low

            elif high_err < reg_err:
                rho_c = rho_c_high

            else:
                diff = diff*2


            rho_c_high = rho_c + diff
            rho_c_low = rho_c - diff
            rho_c_low = max(5000, rho_c_low)

        elif all(err < 0 for err in all_err):
            diff = rho_c_high - rho_c

            if low_err > reg_err:
                rho_c = rho_c_low

            elif high_err > reg_err:
                rho_c = rho_c_high

            else:
                diff = diff*2

            rho_c_high = rho_c + diff
            rho_c_low =  rho_c - diff
            rho_c_low = max(5000, rho_c_low)

        else:
            if  (reg_err * low_err) < 0:
                diff = diff/5
                rho_c_high = rho_c + diff
                rho_c_low = rho_c - diff

        if i > 30:
            break

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
        array2D[index + 1] = star.properties[variable].data(deriv)

    print("Writing star:", name)
    data.array2D2txt(array2D, ["radius"] + save_variable, name)
