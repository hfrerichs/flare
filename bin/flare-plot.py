#!/usr/bin/env PYTHON
import argparse

import matplotlib.pyplot as plt
import flare.plot_data



#===============================================================================
# PLOT command: display data files
#===============================================================================
def plot():
    parser = argparse.ArgumentParser(prog="flare plot")
    parser.add_argument("data_file")
    parser.add_argument("data_key", help="data key as listed by 'flare info'")
    parser.add_argument("-g", "--geometry_file", help="alternative geometry file")
    plot_functions = ['tricontourf', 'scipy_griddata_pcontourf', 'matplotlib_griddata_pcontourf']
    parser.add_argument("-pf", "--plot_function", choices=plot_functions, default=plot_functions[0], help="select plot function for unstructured data")
    parser.add_argument("-D", "--dimension", type=int, help="set plot dimension")
    args   = parser.parse_args()
    flare.plot_data.plot_data(args.data_file, args.data_key, args.geometry_file, plot_function=args.plot_function, dimension=args.dimension)
    plt.show()
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    plot()
