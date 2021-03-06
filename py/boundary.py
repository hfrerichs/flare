from backend import boundary_interface

import matplotlib.pyplot as plt


class Boundary():
    def __init__(self):
        self.n = boundary_interface.num_boundaries()


    def plot_slice(self, phi=0.0):
        for i in range(self.n):
            boundary_interface.get_slice(i+1, phi)
            if boundary_interface.n > 0:
                x = boundary_interface.x
                plt.plot(x[:,0], x[:,1])
