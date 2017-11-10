from copy import copy
import backend

import matplotlib.pyplot as plt



class Separatrix():
    def __init__(self, ix, stop_at_boundary=False):
        backend.separatrix_interface.generate(ix, stop_at_boundary)
        self.b_lc = copy(backend.separatrix_interface.b_lc)
        self.b_rc = copy(backend.separatrix_interface.b_rc)
        self.b_ld = copy(backend.separatrix_interface.b_ld)
        self.b_rd = copy(backend.separatrix_interface.b_rd)

    def plot(self, *args, **kwargs):
        p1, = plt.plot(self.b_lc[:,0], self.b_lc[:,1], *args, **kwargs)
        if 'label' in kwargs: del kwargs['label']
        plt.plot(self.b_rc[:,0], self.b_rc[:,1], p1.get_color(), *args, **kwargs)
        plt.plot(self.b_ld[:,0], self.b_ld[:,1], p1.get_color(), *args, **kwargs)
        plt.plot(self.b_rd[:,0], self.b_rd[:,1], p1.get_color(), *args, **kwargs)
