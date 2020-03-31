"""
The goal of this code is to demonstrate the current status of the star
solver and how to use it. It loops through a star somewhat similar to
the sun and will output a couple of plots to demonstrate.
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import stellar_properties as star
import Use_Data as data


def example_star(starname=""):
    """
    Creates an example star.
    the set of DE's. Afterwards it creates plots of some of the
    desired values to show that it works as expected
    """



    data = pd.read_csv("lowmass_star.txt")



    radius = 6e5*1000
    steps = 7000
    step_size = 1400
    max_step = 1000000

    steps_sizes = np.array([])
    error_list = np.array([])

    #sun_like_star = star.Star(
    #    step_size=step_size,
    #    max_step=max_step,
    #    X=0.60,
    #    Y=0.38,
    #    Z=0.05,
    #    cent_density=58556.,
    #    cent_opticaldepth=0,
    #    cent_temperature=8.23544e6,
    #    cent_radii=700,
    #    error_thresh=1e-5)
    sun_like_star = star.Star(
        step_size=step_size,
        max_step=max_step,
        X=0.60,
        Y=0.38,
        Z=0.05,
        cent_density=99172.323,
        cent_opticaldepth=0,
        cent_temperature=1.8e7,
        cent_radii=700,
        error_thresh=1e-5)

    run = True
    radii_step = 0
    sun_like_star.solve()

    n = 0
    data['rho(r)'] = data['rho(r)'] * (100/1)**3 * (1/1000)
    data['M(r)'] = data['M(r)'] /1000
    data['L(r)'] = data['L(r)'] * 1e-7
    data['r'] = data['r']/100
    data['kappa'] = data['kappa'] * (1/100)**2 * (1000/1)
    data['drho/dr'] = data['drho/dr'] * (100/1)**4 * (1/1000)
    data['dM/dr'] = data['dM/dr'] /1000*(100/1)
    data['dL/dr'] = data['dL/dr'] * 1e-7*(100/1)
    items = ["density", "temperature", "mass", "luminosity"]#, "opacity"]
    comps = ['drho/dr', 'dT/dr', 'dM/dr', 'dL/dr']#, 'kappa']
    comps = ['rho(r)', 'T(r)', 'M(r)', 'L(r)']#, 'kappa']
    # This nixt line is needed for saving data to a text file
    array2D = [[] for i in range(len(items))]
    fig, ax = plt.subplots(2,3, dpi=300)
    n = [0, 0, 0, 1, 1]
    m = [0, 1, 2, 0, 1]
    for item, comp, nn, mm in zip(items,comps, n, m):
        print(nn, mm)
        array2D[0] = sun_like_star.properties[item].data(0)
        print(item, len(array2D[0]))
        ax[nn,mm].plot(
            sun_like_star.properties['radius'],
            array2D[0],  # Plots the value of item
            label=item)
        ax[nn,mm].plot(
            data['r'],
            abs(data[comp]),  # Plots the value of item
            label="true")
        ax[nn,mm].legend()
        #ax[nn,mm].set_yscale('log')
    plt.savefig("compare.png")
    plt.close()
    print(sun_like_star.properties['density'].data(0))


    # Get user decision on saving star data
    save = input("Save this star? (y/n): ")
    # Save the data to a text file
    if save == 'y': data.array2D2txt(array2D, items, starname)


if __name__ == '__main__':
    example_star()
