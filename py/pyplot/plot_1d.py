import os
import numpy             as np
import matplotlib.pyplot as plt


def plot_1d(data_file, *args, **kwargs):
    # check if data file exists
    if not os.path.isfile(data_file):
        print "error: data file {} not found!".format(data_file)
        return

    # load data file and check dimensions
    d = np.loadtxt(data_file, dtype='float')
    if d.ndim != 2:
        print "error: unexpected dimension of data array!"
        return

    fill = kwargs.get('fill', False)

    if fill:
        plt.fill(d[:,0], d[:,1], *args, **kwargs)
    else:
        plt.plot(d[:,0], d[:,1], *args, **kwargs)
