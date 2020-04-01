import argparse as arg
from make_star import make_star

def main(args):
    file = open(args.fileName, 'r')
    file_lines = file.readlines()
    file_lines = [file for file in file_lines if '#' not in file]
    for line in file_lines:
        line = line.replace("\n","").split(", ")
        print(line)
        name = "Tc_{:.2e}_rhoc_guess{:.2e}_Core_{}".format(
                float(line[0]), float(line[1]), line[2], line[3])

        make_star(line[0], line[1], line[2], name)

if __name__ == '__main__':
    parser = arg.ArgumentParser(description = "Plots Stars!")
    parser.add_argument('fileName',
                        help='Enter the file name that contains the list of stars that you want to run')
    args = parser.parse_args()

    main(args)
