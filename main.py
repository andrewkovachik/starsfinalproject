import argparse as arg
from multiprocessing import Pool
from make_star import make_star

def unpack(line):
    print(line)
    line = line.replace("\n","").split(", ")
    name = "Tc_{:.2e}_rhoc_guess{:.2e}_Core_{}_Type_{}".format(
            float(line[0]), float(line[1]), line[2], line[3])

    args = (float(line[0]), float(line[1]), line[2], name)
    return make_star(*args)

def main(args):
    pool = Pool()
    file = open(args.fileName, 'r')
    file_lines = file.readlines()
    file_lines = [file for file in file_lines if '#' not in file]
    last_rho_c = 0
    if args.parallel:
        print("Running Parallel")
        results = pool.map(unpack, file_lines)

    else:
        for line in file_lines:
            line = line.replace("\n","").split(", ")

            print()
            print(line)
            if last_rho_c and args.adaptive:
                name = "Tc_{:.2e}_rhoc_guess{:.2e}_Core_{}_Type_{}".format(
                        float(line[0]), last_rho_c, line[2], line[3])
            else:
                name = "Tc_{:.2e}_rhoc_guess{:.2e}_Core_{}_Type_{}".format(
                        float(line[0]), float(line[1]), line[2], line[3])

            try:
                if last_rho_c and args.adaptive:
                    last_rho_c = make_star(float(line[0]), last_rho_c, line[2], name)
                else:
                    last_rho_c = make_star(float(line[0]), float(line[1]), line[2], name)
            except:
                print("Failed making star %s"%(name))


if __name__ == '__main__':
    parser = arg.ArgumentParser(description = "Plots Stars!")
    parser.add_argument('fileName',
                        help='Enter the file name that contains the list of stars that you want to run')
    parser.add_argument('--parallel',
                        action='store_true',
                        help='Run multiple stars at the same time')
    parser.add_argument('--adaptive',
                        action='store_true',
                        help='Use the previous solutions rho_c as the guess for this one. May speed up if stars change linearly')
    args = parser.parse_args()

    main(args)
