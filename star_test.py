"""
The goal of this code is to demonstrate the current status of the star
solver and how to use it. It loops through a star somewhat similar to
the sun and will output a couple of plots to demonstrate.
"""
import matplotlib.pyplot as plt
import numpy as np
import stellar_properties as star

def example_star():
    """
    Creates an example star, prints out some of the steps as it solves
    the set of DE's. Afterwards it creates plots of some of the
    desired values to show that it works as expected
    """

    expected_radii = 300000000
    steps = 5000
    step_size = expected_radii / steps

    sun_like_star = star.Star(
        step_size=step_size,
        X=0.55,
        Y=0.43,
        Z=0.02,
        cent_density=162200,
        cent_opticaldepth=2/3,
        cent_temperature=1.5*10**7,
        cent_radii=0.01
        )

    radii_steps = np.arange(steps)
    for radii_step in radii_steps:
        sun_like_star.step_de()
        sun_like_star.step_non_de()

        if radii_step % 200 == 0:
            print(sun_like_star)

    items = ["density", "temperature", "mass", "luminosity", "opticaldepth"]
    for item in items:
        plt.plot(
            radii_steps * step_size / 696342000, # Units of solar radius
            sun_like_star.properties[item].val[0, :-1], # Plots the value of item
            label=item)
        plt.xlabel("Solar Radii")
        plt.title(item)
        plt.savefig("%s.png" % (item))
        plt.close()

if __name__ == '__main__':
    example_star()
