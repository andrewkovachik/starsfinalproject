import os
from pathlib import Path
import Use_Data as data

def plotdata(toPlot, filename, folder="Main_Sequences_txt_Files"):
	"""
	plotdata takes in an array with the names of columns that are 
	to be plotted, a filename with the data, and optionally a folder
	which contains the files

	Args:
	    toPlot (np.array): array of names of values to be plotted
	    filename (str): identifier for file
	    folder (str): folder name without forward slash
	"""

	# Create path to file
	filepath = folder + "/" + filename
	arr = data.txt2array2D("/Users/updownleftright/Documents/GitHub/starsfinalproject/Main_Sequences_txt_Files/star")
	cols = [0]
	ln = len(toPlot)

	for i in range(ln):
		for j in range(ln):
			test = 1
		if j == ln - 1: 
			print("j is ln ", i)


plotdata([1,2,3,4,5],"star.txt")


