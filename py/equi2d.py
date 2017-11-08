import sys
import numpy as np
import matplotlib.pyplot as plt

from backend  import equi2d_interface
from boundary import *


LOAD_ERROR = ['',
    'invalid format string given!',
    'guessing equilibrium format failed!',
    'loading data file failed!']


class LoadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def load(filename, format):
    ierr = equi2d_interface.load(filename, format)
    if ierr > 0:
        raise(LoadError(LOAD_ERROR[ierr]))


def view(filename, format):
    try:
        load(filename, format)
    except LoadError as e:
        print 'error: ', e.value
        print 'use "-f FORMAT" to provide equilibrium format'
        sys.exit(2)


    # sample poloidal flux
    rmin, rmax, zmin, zmax = equi2d_interface.get_domain()
    nr  = 128
    nz  = 128
    r   = np.linspace(rmin, rmax, nr)
    z   = np.linspace(zmin, zmax, nz)
    psi = equi2d_interface.get_psi(r, z)


    # plot poloidal flux
    levels = np.linspace(np.min(psi), np.max(psi), 64)
    plt.contourf(r, z, psi, levels=levels)
    plt.xlabel("Major radius [cm]")
    plt.ylabel("Vertical position [cm]")


    # plot boundary (if available)
    B = Boundary()
    B.plot_slice(0.0)


    # add color bar
    cbar = plt.colorbar()
    cbar.set_label("Poloidal flux [weber]")
    plt.show()


def init(filename, format):
    try:
        load(filename, format)
    except LoadError as e:
        print 'error: ', e.value
        print 'use "-f FORMAT" to provide equilibrium format'
        sys.exit(2)


    ierr = equi2d_interface.init(filename, format)
