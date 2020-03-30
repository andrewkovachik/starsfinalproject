import Use_Data as data
import matplotlib.pyplot as plt
from pathlib import Path

def plotdata(toPlot, filename, folder="Main_Sequences_txt_Files", 
	title="", xtitle="", ytitle="", xlim=[0,0], ylim=[0,0]):
	"""
	plotdata takes in an array with the names of columns that are 
	to be plotted, a filename with the data, and optionally a folder
	which contains the files as well as many plotting options


	Args:
	    toPlot (np.array): array of names of values to be plotted
	    filename (str): identifier for file
	    folder (str): folder name without forward slash
	"""

	# Get data and header from file
	arr, header = data.txt2array2D(filepath)

	# Set up for the plotting loop
	plotlen = len(toPlot)
	headlen = len(header)
	exists = False
	label = -1

	# Plotting loop, starts by checking if toPlot values are valid
	for i in range(plotlen):
		for j in range(headlen):
			if toPlot[i] == header[j]:
				exists = True
				label = j
				break

		if not exists: # If toPlot value doesn't exist tell the user
			print("There's no", toPlot[i], "in the", filename, "file")
		
		else: # Plot data
			plt.plot(arr[0], arr[i], label=header[label])
		
		# Reset to find next toPlot label
		exists = False
		label = -1

	# Optional plotting things
	if title  != "": plt.title(title)
	if xtitle != "": plt.xlabel(xtitle)
	if ytitle != "": plt.ylabel(ytitle)
	if xlim != [0,0]: plt.xlim(xlim)
	if ylim != [0,0]: plt.ylim(ylim)
	plt.legend()
	plt.show()


plotdata(["density", "temperature", "mass", "luminosity", "opticaldepth"],"star")
