from flare import *

import numpy             as np
import matplotlib.pyplot as plt


def plot_lc(grid_file, data_file, qplot='lc'):
    # load grid file
    G = Grid(grid_file)
    if G.layout != STRUCTURED:
        print "error: cannot visualize unstructured grid!"
        return

    # load data file
    d = np.loadtxt(data_file, dtype='float')
    if d.ndim != 2:
        print "error: unexpected dimension of data array!"
        return

    # backward connection length [m]
    lc_bwd = abs(d[:,0].reshape(G.n2, G.n1)) / 100.0
    # forward connection length [m]
    lc_fwd = abs(d[:,1].reshape(G.n2, G.n1)) / 100.0
    # total connection length [m]
    lc     = lc_bwd + lc_fwd


    # select quantity for plotting
    q = {
        'lc':     lc,
        'lcs':    np.minimum(lc_bwd, lc_fwd),
        'lc_bwd': lc_bwd,
        'lc_fwd': lc_fwd
    }.get(qplot)
    if q is None:
        print "error: invalid quantity '{}'!".format(qplot)
        return
    qmin = np.nanmin(q)
    qmax = np.nanmax(q)


    # plot options
    levels = np.linspace(qmin, qmax, 64)
    plt.contourf(G.x1, G.x2, q, vmin=qmin, vmax=qmax, levels=levels)
    plt.xlabel(G.x1_label)
    plt.ylabel(G.x2_label)
    cbar = plt.colorbar()
    cbar.set_label('Connection Length [m]')
    plt.show()
