from flare import *

import numpy             as np
import matplotlib.pyplot as plt


def plot_2d(grid_file, data_file, icol=0):
    # load grid file
    G = Grid(grid_file)
    if G.layout != STRUCTURED:
        print "error: cannot visualize unstructured grid!"
        return

    # load data file
    d = np.loadtxt(data_file, dtype='float')
    d = d.reshape(d.size,1)
    if d.ndim != 2:
        print "error: unexpected dimension of data array!"
        print d.ndim
        return

    # select quantity for plotting
    q = d[:,icol].reshape(G.n2, G.n1)
    qmin = np.nanmin(q)
    qmax = np.nanmax(q)


    # plot options
    levels = np.linspace(qmin, qmax, 64)
    plt.contourf(G.x1, G.x2, q, vmin=qmin, vmax=qmax, levels=levels)
    plt.xlabel(G.x1_label)
    plt.ylabel(G.x2_label)
    cbar = plt.colorbar()
    plt.show()
