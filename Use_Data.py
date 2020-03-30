import os
from pathlib import Path


def array2D2txt(array, header=[], filename="", folder="Star_Files"):
    """
    results2txt takes in a 2D array, a filename without an extension, and a
    folder path to create a .txt file with the values of the 2D array in
    the folder with the given filename, and if no folder or filename are
    specified then a default value is given for both

    Args:
        array (np.array): data to be written to a file
        filename (str): identifier for file
        folder (str): folder name without forward slash
    """
    # Checks if there is a folder and makes one if there is not
    Path(folder).mkdir(parents=True, exist_ok=True)

    # Create path to file
    filepath = folder + "/" + filename
    # Make defualt name
    default = "star"

    # If no file name is set come name it the default
    if filepath == folder + "/":
        filepath += default

    # Conditions for if a file name is given
    if os.path.exists(filepath):
        # If filename exists rename it with a number appended
        counter = 0
        filepath += "_{}.txt"

        while os.path.exists(filepath.format(counter)):
            # Loop until a name exists without some number
            counter +=1
        filepath = filepath.format(counter)

    # Open a file with filename and write in the values that are tab seperated
    text_file = open(filepath, "w")

    # Write header if there is one
    if header != []:
        ln = len(header)
        for i in range(ln):
            n = text_file.write(str(header[i]) + "\t")
        n = text_file.write("\n")

    # Write content
    for i in range(len(array[0])):
        for j in array:
            n = text_file.write(str(j[i]) + "\t")
        n = text_file.write("\n")

    # Close text file
    text_file.close()


def txt2array2D(filepath):
    """
    txt2array2D takes in a file path and outputs a 2D array with the vaules
    from the txt file given that each line is tab separated

    Args:
        filepath (str): file path to search for and read in

    Return:
        (array): Text file data in the form of a 2D array
    """
    # Open the text file
    text_file = open(filepath, "r")

    # Get the number of rows and columns from text file and make a 2D array
    file = text_file.read().split("\n")
    cols = len(file[0].split("\t")) - 1
    rows = len(file) - 1
    array = [[0.0 for i in range(rows)] for j in range(cols)]

    # Gets header name line from file
    if isinstance(file[0], str):
        header = file[0].split("\t")
        header.pop(len(header) - 1)
        file.pop(0)
        rows = rows - 1

    # Put the values from the text file into the 2D array
    for i in range(rows):
        line = file[i].split("\t")
        for j in range(cols):
            array[j][i] = float(line[j])

    text_file.close()

    return array, header
