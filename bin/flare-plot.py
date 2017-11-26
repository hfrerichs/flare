#!/usr/bin/env python
import argparse

import matplotlib.pyplot as plt
import flare.pyplot      as fplt



#===============================================================================
# PLOT command: display data files
#===============================================================================
def plot():
    parser = argparse.ArgumentParser(prog="flare plot")
    parser.add_argument("data_key", help="data key as listed by 'flare info'")
    parser.add_argument("-g", "--geometry_file", help="alternative geometry file")
    plot_functions = ['tricontourf', 'scipy_griddata_pcontourf', 'matplotlib_griddata_pcontourf']
    parser.add_argument("-pf", "--plot_function", choices=plot_functions, default=plot_functions[0], help="select plot function for unstructured data")
    parser.add_argument("data_file")
    args   = parser.parse_args()
    fplt.plot_data(args.data_file, args.data_key, args.geometry_file, plot_function=args.plot_function)
    plt.show()
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    plot()
