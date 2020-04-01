import Use_Data as data
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import rc
import argparse as arg

#rc('text', usetex=True) # For latex in plots


def plotdata(toPlot,
             filename,
             folder="Star_Files",
             save=False,
             title="",
             xtitle="",
             ytitle="",
             xlim=[0, 0],
             ylim=[0, 0]):
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
    exists = False
    label = -1

    for item in toPlot:
        if item in header:
            plt.plot(
                arr[0,:-1],
                arr[header.index(item),:-1] / max(arr[header.index(item)]),
                label=item)

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


if __name__ == '__main__':
    parser = arg.ArgumentParser(description="Plots Stars!")
    parser.add_argument(
        'fileName',
        help='Enter the file name that contains the star you want to plot',
        default='star')

    args = parser.parse_args()
    toPlot = ['density', 'temperature', 'opticaldepth', 'mass', 'luminosity']
    filename = args.fileName
    title = "Test Plot"
    xtitle = "Test X Axis"
    ytitle = "Test Y Axis"
    plotdata(
        toPlot, filename, save=True, title=title, xtitle=xtitle, ytitle=ytitle)
