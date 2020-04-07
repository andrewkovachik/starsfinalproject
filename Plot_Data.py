import argparse as arg
import numpy as np
import Use_Data as data
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import rc
from os import listdir

#rc('text', usetex=True) # For latex in plots


def plotdata(toPlot,
             filename,
             folder="Star_Files",
             save=False,
             title="",
             xtitle="",
             ytitle="",
             xlim=[0, 0],
             ylim=[0, 0],
             divmax=True):
    """
	plotdata takes in an array with the names of columns that are 
	to be plotted, a filename with the data, and optionally a folder
	which contains the files as well as many plotting options

	Args:
			toPlot (np.array): array of names of values to be plotted
			filename (str): identifier for file
			folder (str): folder name without forward slash
	"""

    filepath = folder + "/" + filename

    # Get data and header from file
    arr, header = data.txt2array2D(filepath)

    # Set up for the plotting loop
    plotlen = len(toPlot)
    headlen = len(header)

    for item in toPlot:
        if item in header:
            pltarr = arr[header.index(item)]
            pltmax = max(pltarr)
            if divmax:
                pltarr = [i / pltmax for i in pltarr]
                plt.plot(arr[0], pltarr, label=item)
            else:
                plt.plot(arr[0], plrarr, label=item)

        else:  # If toPlot value doesn't exist tell the user
            print("There's no", item, "in the", filename, "file")

    # Optional plotting things
    if title != "": plt.title(title)
    if xtitle != "": plt.xlabel(xtitle)
    if ytitle != "": plt.ylabel(ytitle)
    if xlim != [0, 0]: plt.xlim(xlim)
    if ylim != [0, 0]: plt.ylim(ylim)
    plt.figure(1, dpi=300)
    plt.legend()
    if save: plt.savefig(filepath[:-4] + ".png")
    if not save: plt.show()
    plt.close()


def plotall(toPlot,
            folder="Star_Files",
            save=False,
            title="",
            xtitle="",
            ytitle="",
            xlim=[0, 0],
            ylim=[0, 0],
            divmax=True):
    """
	plotdata takes in an array with the names of columns that are 
	to be plotted and optionally a folder which contains the files 
	as well as many plotting options for all the files in the folder 
	that end in .txt
	
	Args:
		toPlot (np.array): array of names of values to be plotted
		folder (str): folder name without forward slash
	"""
    files = listdir(folder + "/")  # Get list of files in folder

    for file in files:
        if file[-4:] == ".txt":  # Only continue for text files
            name = file.split("_")  # Get file name to make title
            title = (name[5] + ' Core, ' + name[7][:2] + ' Star, T$_c$ = ' +
                     name[1] + ', $\\rho_c$ = ' + name[3][5:])

            # Plot data with the proper title
            plotdata(
                toPlot,
                file,
                save=True,
                title=title,
                xtitle='Radius (m)',
                ytitle='Components / Max Value')

            print("Plotted star:", title)
            print()


def plotmass_lum(folder="Star_Files"):
    files = listdir(folder + "/")  # Get list of files in folder
    # Make temp and lum arrays
    luminosityH = []
    luminosityHe = []
    luminosityC = []
    radiusH = []
    radiusHe = []
    radiusC = []

    for file in files:  # Loop through the files
        if file[-4:] == ".txt":  # If it's a text file
            # Get the text file data
            arr, header = data.txt2array2D(folder + "/" + file)

            name = file.split("_")

            if name[5] == "Hydrogen":
                luminosityH.append(arr[header.index("luminosity")][-1])
                radiusH.append(arr[header.index("radius")][-1])

            if name[5] == "Helium":
                luminosityHe.append(arr[header.index("luminosity")][-1])
                radiusHe.append(arr[header.index("radius")][-1])

            if name[5] == "Carbon":
                luminosityC.append(arr[header.index("luminosity")][-1])
                radiusC.append(arr[header.index("radius")][-1])

            print("Got star", file)

    luminosityH = np.array(luminosityH)/3.828e26
    luminosityHe = np.array(luminosityHe)/3.828e26
    luminosityC = np.array(luminosityC)/3.828e26
    radiusH = np.array(radiusH)/6.957e8
    radiusHe = np.array(radiusHe)/6.957e8
    radiusC = np.array(radiusC)/6.957e8

    plt.scatter(
        np.log10(radiusH),
        np.log10(luminosityH),
        c="Blue",
        marker="^",
        label="Hydrogen",
        alpha=0.7)
    plt.scatter(
        np.log10(radiusHe),
        np.log10(luminosityHe),
        c="Red",
        marker="s",
        label="Helium",
        alpha=0.7)
    plt.scatter(
        np.log10(radiusC),
        np.log10(luminosityC),
        c="Green",
        label="Carbon",
        alpha=0.7)
    plt.title("Main Sequence")
    plt.xlabel("log(Luminosity (L$_\\odot$))")
    plt.ylabel("log(Radius (R$_\\odot$))")
    plt.legend()
    plt.show()
def plotmass_radius(folder="Star_Files"):
    files = listdir(folder + "/")  # Get list of files in folder
    # Make temp and lum arrays
    massH = []
    massHe = []
    massC = []
    radiusH = []
    radiusHe = []
    radiusC = []

    for file in files:  # Loop through the files
        if file[-4:] == ".txt":  # If it's a text file
            # Get the text file data
            arr, header = data.txt2array2D(folder + "/" + file)

            name = file.split("_")

            if name[5] == "Hydrogen":
                massH.append(arr[header.index("mass")][-1])
                radiusH.append(arr[header.index("radius")][-1])

            if name[5] == "Helium":
                massHe.append(arr[header.index("mass")][-1])
                radiusHe.append(arr[header.index("radius")][-1])

            if name[5] == "Carbon":
                massC.append(arr[header.index("mass")][-1])
                radiusC.append(arr[header.index("radius")][-1])

            print("Got star", file)

    massH = np.array(massH)/1988e30
    massHe = np.array(massHe)/1988e30
    massC = np.array(massC)/1988e30
    radiusH = np.array(radiusH)/6.957e8
    radiusHe = np.array(radiusHe)/6.957e8
    radiusC = np.array(radiusC)/6.957e8

    plt.scatter(
        np.log10(radiusH),
        np.log10(massH),
        c="Blue",
        marker="^",
        label="Hydrogen",
        alpha=0.7)
    plt.scatter(
        np.log10(radiusHe),
        np.log10(massHe),
        c="Red",
        marker="s",
        label="Helium",
        alpha=0.7)
    plt.scatter(
        np.log10(radiusC),
        np.log10(massC),
        c="Green",
        label="Carbon",
        alpha=0.7)
    plt.title("Main Sequence")
    plt.xlabel("log(Mass (M$_\\odot$))")
    plt.ylabel("log(Radius (R$_\\odot$))")
    plt.legend()
    plt.show()

def plotmain(folder="Star_Files"):
    files = listdir(folder + "/")  # Get list of files in folder
    # Make temp and lum arrays
    temperatureH = []
    luminosityH = []
    temperatureHe = []
    luminosityHe = []
    temperatureC = []
    luminosityC = []

    for file in files:  # Loop through the files
        if file[-4:] == ".txt":  # If it's a text file
            # Get the text file data
            arr, header = data.txt2array2D(folder + "/" + file)

            name = file.split("_")

            if name[5] == "Hydrogen":
                for item in header:  # Find the temp and lum columns
                    if item == 'temperature':  # Add to temp
                        temperatureH.append(arr[header.index(item)][-1])

                    elif item == 'luminosity':  # Add to lum
                        luminosityH.append(
                            arr[header.index(item)][-1] / 3.828e26)

            if name[5] == "Helium":
                for item in header:  # Find the temp and lum columns
                    if item == 'temperature':  # Add to temp
                        temperatureHe.append(arr[header.index(item)][-1])

                    elif item == 'luminosity':  # Add to lum
                        luminosityHe.append(
                            arr[header.index(item)][-1] / 3.828e26)

            if name[5] == "Carbon":
                for item in header:  # Find the temp and lum columns
                    if item == 'temperature':  # Add to temp
                        temperatureC.append(arr[header.index(item)][-1])

                    elif item == 'luminosity':  # Add to lum
                        luminosityC.append(
                            arr[header.index(item)][-1] / 3.828e26)

            print("Got star", file)

    plt.scatter(
        np.log10(temperatureH),
        np.log10(luminosityH),
        c="Blue",
        marker="^",
        label="Hydrogen",
        alpha=0.7)
    plt.scatter(
        np.log10(temperatureHe),
        np.log10(luminosityHe),
        c="Red",
        marker="s",
        label="Helium",
        alpha=0.7)
    plt.scatter(
        np.log10(temperatureC),
        np.log10(luminosityC),
        c="Green",
        label="Carbon",
        alpha=0.7)
    plt.title("Main Sequence")
    plt.ylabel("log(Luminosity (L$_\\odot$))")
    plt.xlabel("log(Surface Temperature (K))")
    plt.legend()
    plt.gca().invert_xaxis()
    ax = plt.gca()

    ax.axvspan(np.log10(35000), np.log10(30000), alpha=0.4, color="xkcd:sky")
    ax.axvspan(
        np.log10(30000), np.log10(10000), alpha=0.4, color="xkcd:powder blue")
    ax.axvspan(
        np.log10(10000), np.log10(7500), alpha=0.4, color="xkcd:pale blue")
    ax.axvspan(
        np.log10(6000), np.log10(7500), alpha=0.4, color="xkcd:off white")
    ax.axvspan(
        np.log10(5200), np.log10(6000), alpha=0.4, color="xkcd:light yellow")
    ax.axvspan(
        np.log10(3700), np.log10(5200), alpha=0.4, color="xkcd:yellow orange")
    ax.axvspan(
        np.log10(1000), np.log10(3700), alpha=0.4, color="xkcd:orange red")

    plt.text(np.log10(12500), -5, '$B$', color='xkcd:dark grey', zorder=100)
    plt.text(
        np.log10((10000 - 7500) / 2 + 7500),
        -5,
        '$A$',
        color='xkcd:dark grey',
        zorder=100)
    plt.text(
        np.log10((7500 - 6000) / 2 + 6000),
        -5,
        '$F$',
        color='xkcd:dark grey',
        zorder=100)
    plt.text(
        np.log10((6000 - 5200) / 2 + 5300),
        -5,
        '$G$',
        color='xkcd:dark grey',
        zorder=100)
    plt.text(
        np.log10((5200 - 3700) / 2 + 3850),
        -5,
        '$K$',
        color='xkcd:dark grey',
        zorder=100)
    plt.text(
        np.log10((3700 - 1000) / 2 + 1000),
        -5,
        '$M$',
        color='xkcd:dark grey',
        zorder=100)

    plt.show()


if __name__ == '__main__':
    parser = arg.ArgumentParser(description="Plots Stars!")
    parser.add_argument(
        'fileName',
        help='Enter the file name that contains the star that you want to plot',
        default="")

    args = parser.parse_args()
    toPlot = ['density', 'temperature', 'opticaldepth', 'mass', 'luminosity']
    #filename = "Tc_1.54e+07_rhoc_guess7.00e+05_Core_Hydrogen_Type_MS.txt"
    title = "Test Plot"
    xtitle = "Test X Axis"
    ytitle = "Test Y Axis"

    if args.fileName == "all":
        plotall(toPlot)
        print("All star plot png's created and saved")

    elif args.fileName == "main":
        plotmain()
        print("Main sequence plotted (No save)")

    elif args.fileName == "lum":
        plotmass_lum()
        print("Plotting Stars as a function of mass and radii")

    elif args.fileName == "mass":
        plotmass_radius()
        print("Plotting Stars as a function of mass and radii")

    else:
        plotdata(
            toPlot,
            args.fileName,
            save=True,
            title=title,
            xtitle=xtitle,
            ytitle=ytitle)
        print("Star:", args.fileName, "plot saved")
    print()
