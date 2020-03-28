"""
The goal of this code is to demonstrate the current status of the star
solver and how to use it. It loops through a star somewhat similar to
the sun and will output a couple of plots to demonstrate.
"""
import matplotlib.pyplot as plt
import numpy as np
import stellar_properties as star
import Use_Data as data


def example_star(starname=""):
    """
    Creates an example star.
    the set of DE's. Afterwards it creates plots of some of the
    desired values to show that it works as expected
    """

    radius = 6e5*1000
    steps = 400
    step_size = 100
    max_step = 1000000

    steps_sizes = np.array([])
    error_list = np.array([])

    sun_like_star = star.Star(
        step_size=step_size,
        max_step=max_step,
        X=0.55,
        Y=0.43,
        Z=0.02,
        cent_density=58556.,
        cent_opticaldepth=0,
        cent_temperature=8.23544e6,
        cent_radii=100,
        error_thresh=10e-5)

    radii_steps = np.arange(steps)
    for radii_step in radii_steps:
        sun_like_star.step_de()
        sun_like_star.step_non_de()

        radii_steps = np.append(radii_steps, sun_like_star.step_size)
        error_list = np.append(error_list, max(sun_like_star.error))


        if radii_step%20 == 0:
            print(sun_like_star)

    n = 0
    items = ["density", "temperature", "mass", "luminosity", "opticaldepth"]
    # This nixt line is needed for saving data to a text file
    array2D = [[] for i in range(len(items))]
    for item in items:
        array2D[n] = sun_like_star.properties[item].val[0, :] / max(
            sun_like_star.properties[item].val[0, :])
        plt.plot(
            sun_like_star.properties['radius'] /
            696342000,  # Units of solar radius
            array2D[n],  # Plots the value of item
            label=item)
    plt.xlabel("Solar Radii")
    plt.legend()
    plt.savefig("RelativeAmps.png")
    plt.close()

    plt.plot(np.arange(len(radii_steps)), radii_steps)
    plt.savefig("RadiusSteps.png")
    plt.close()
    plt.plot(np.arange(len(error_list)), error_list)
    plt.savefig("Error.png")
    plt.close()
    n += 1

    # Get user decision on saving star data
    save = input("Save this star? (y/n): ")
    # Save the data to a text file
    if save == 'y': data.array2D2txt(array2D, items, starname)


if __name__ == '__main__':
    example_star()
