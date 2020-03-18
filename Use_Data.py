import os
from pathlib import Path


def array2D2txt(array, filename="", folder="Main_Sequences_txt_Files"):
    """
    results2txt takes in a 2D array, a filename without an extension, and a
    folder path to create a .txt file with the values of the 2D array in
    the folder witht he given filename, and if no folder or filename are
    specified then a default value for both are given

    Args:
        array (np.array): data to be written to a file
        filename (str): identifier for file
        folder (str): folder name without forward slash
    """
    # Checks if there is a folder and makes one if there is not
    Path(folder).mkdir(parents=True, exist_ok=True)

    # Create path to file
    filepath = folder + "/" + filename

    # If no file name is set come up with a standard
    if filepath == folder + "/":
        # gets the number of files in the folder for naming purposes
        num_files = len(os.listdir(folder))
        # file is renamed based on the number of files in the folder
        filepath += "star_" + str(num_files) + ".txt"

    # Conditions for if a file name is given
    elif os.path.exists(filepath):
        # If filename exists rename it with
        filepath.replace(".txt", "_1.txt")

    # Open a file with filename, write in the values that are tab seperated,
    #  and then close the file
    text_file = open(filepath, "w")
    for i in range(len(array[0])):
        for j in array:
            n = text_file.write(str(j[i]) + "\t")
        n = text_file.write("\n")
    text_file.close()


def txt2array2D(filepath):
    """
    txt2array2D takes in a file path and outputs a 2D array with the vaules
    from the txt file given that each line is tab separated

    Args:
        filepath (str): file path to search for and read in

    Return:
        (array): Text file data in the form of a 2d array
    """
    # Open the text file
    text_file = open(filepath, "r")

    # Get the number of rows and columns from text file and make a 2D array
    file = text_file.read().split("\n")
    cols = len(file[0].split("\t")) - 1
    rows = len(file) - 1
    array = [[0.0 for i in range(rows)] for j in range(cols)]

    # Put the values from the text file into the 2D array
    for i in range(rows):
        line = file[i].split("\t")
        for j in range(cols):
            array[j][i] = float(line[j])

    text_file.close()

    return array
