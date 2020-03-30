"""
The goal of this code is to demonstrate the current status of the star
solver and how to use it. It loops through a star somewhat similar to
the sun and will output a couple of plots to demonstrate.
"""
import matplotlib.pyplot as plt
import numpy as np
import stellar_properties as star
import Use_Data as data

#Constants
G = 6.6741 * 10**-11
Mp = 1.67 * 10**-27
Kb = 1.381 * 10**-23

#Star mass and radius inital guesses
M=1.989*10**30
R=696340000

            
def example_star(starname=""):
    """
    Creates an example star, prints out some of the steps as it solves
    the set of DE's. Afterwards it creates plots of some of the
    desired values to show that it works as expected
    """

    expected_radii = R
    steps = 100
    step_size = expected_radii / steps

    sun_like_star = star.Star(
        step_size=step_size,
        X=0.55,
        Y=0.43,
        Z=0.02,
        cent_density=162200,
        cent_opticaldepth=2/3,
        #Initial guess for total mass and radius of star based on desired position on HR diagram
        #Define M, R above!
        cent_temperature=2*G*M*Mp/(3*Kb*R),
        cent_radii=0.01
        )

    radii_steps = np.arange(steps)
    for radii_step in radii_steps:
        sun_like_star.step_de()
        sun_like_star.step_non_de()

        print(sun_like_star)

    items = ["density", "temperature", "mass", "luminosity", "opticaldepth"]
    # This nixt line is needed for saving data to a text file
    array2D = [[] for i in range(len(items) + 1)]
    array2d[0] = radii_steps * step_size / 696342000, # Units of solar radius
    n = 1

    for item in items:
        array2D[n] = sun_like_star.properties[item].val[0, :-1]
        plt.plot(
            array2D[0], # Units of solar radius
            array2D[n], # Plots the value of item
            label=item)
        plt.xlabel("Solar Radii")
        plt.title(item)
        #plt.savefig("%s.png" % (item)) # Don't need because of plotting script
        plt.close()
        n += 1

    # Get user decision on saving star data
    save = input("Save this star? (y/n): ")
    # Save the data to a text file
    if save == 'y': data.array2D2txt(array2D, [radius] + items, starname)

if __name__ == '__main__':
    example_star()
