import os
import numpy             as np
import matplotlib.pyplot as plt


def plot_1d(data_file, *args, **kwargs):
    if not os.path.isfile(data_file): return
    # load data file
    d = np.loadtxt(data_file, dtype='float')
    if d.ndim != 2:
        print "error: unexpected dimension of data array!"
        return

    plt.plot(d[:,0], d[:,1], *args, **kwargs)
