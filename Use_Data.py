import os
from pathlib import Path

# results2txt takes in a 2D array, a filename without an extension, and a 
#  folder path to create a .txt file with the values of the 2D array in
#  the folder witht he given filename, and if no folder or filename are
#  specified then a default value for both are given
def array2D2txt(array, filename="", folder="Main_Sequences_txt_Files"):
	# Checks if there is a folder and makes one if there is not
	Path(folder).mkdir(parents=True, exist_ok=True)
	
	# Create path to file
	filename = folder + "/" + filename

	# If no file name is set come up with a standard
	if filename == folder + "/":
		# gets the number of files in the folder for naming purposes
		num_files = len(os.listdir(folder))
		# file is renamed based on the number of files in the folder
		filename += "star_" + str(num_files) + ".txt"

	# Conditions for if a file name is given
	else:
		# If filename exists rename it with
		if os.path.exists(path): filename += "_1.txt"

	# Open a file with filename, write in the values that are tab seperated,
	#  and then close the file
	text_file = open(filename, "w")
	for i in range(len(array[0])):
		for j in array:
			n = text_file.write(str(j[i]) + "\t")
		n = text_file.write("\n")
	text_file.close()

# txt2array2D takes in a file path and outputs a 2D array with the vaules
#  from the txt file given that each line is tab separated
def txt2array2D(filepath):
	# Open the text file
	text_file = open(filepath,"r") 

	# Get the number of rows and columns from text file and make a 2D array
	file = text_file.read().split("\n")
	cols = len(file[0].split("\t")) - 1
	rows = len(file) - 1
	array = [[0.0 for i in range(rows) ] for j in range(cols) ]

	# Put the values from the text file into the 2D array
	for i in range(rows):
		line = file[i].split("\t")
		for j in range(cols):
			array[j][i] = float(line[j])

	text_file.close()

	return array
