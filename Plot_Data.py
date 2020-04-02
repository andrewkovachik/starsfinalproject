import argparse as arg
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
    if save: plt.savefig(filepath + ".png")
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
            print(name)
            title = (name[5] + ' Core, T_c = ' + name[1] + ', rho_c = ' + name[3][5:])

            # Plot data with the proper title
            plotdata(
                toPlot,
                file,
                save=True,
                title=title,
                xtitle='Radius (m)',
                ytitle='Components / Max Value')


if __name__ == '__main__':
    parser = arg.ArgumentParser(description="Plots Stars!")
    parser.add_argument(
        '-fileName',
        help='Enter the file name that contains the star that you want to plot',
        default=""
    )
    args = parser.parse_args()

    toPlot = ['density', 'temperature', 'opticaldepth', 'mass', 'luminosity']
    #filename = "Tc_1.54e+07_rhoc_guess7.00e+05_Core_Hydrogen_Type_MS.txt"
    title = "Test Plot"
    xtitle = "Test X Axis"
    ytitle = "Test Y Axis"

    if args.fileName == "":
        plotall(toPlot)
    else:
        plotdata(
            toPlot,
            args.fileName,
            save=True,
            title=title,
            xtitle=xtitle,
            ytitle=ytitle)
