import matplotlib.pyplot as plt

from PySide.QtGui import *
from PySide.QtCore import *

import flare
from main_ui import Ui_MainWindow

import numpy as np

class FlareUi(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(FlareUi, self).__init__()
        self.setupUi(self)
        self.assignWidgets()
        self.show()

        flare.backend.init()
        flare.backend.init_bfield()

    def assignWidgets(self):
        self.run.clicked.connect(self.run_proc)

    def run_proc(self):
        a = [110.0, 0.0, 0.0]
        #flare.backend.trace_bline_py(a, 1000)

        n  = 64
        x0 = np.zeros([n,n,3])
        for i in range(n):
            for j in range(n):
                x0[i,j,0] = 100.0 + 60.0*i/(n-1)
                x0[i,j,1] =-140.0 + 60.0*j/(n-1)
                x0[i,j,2] =   0.0


        m = n**2
        x0_list = x0.reshape(m,3)
        d = flare.backend.connection_length_py(x0_list, m)

#        for i in range(16):
#            print x0[i,:], d[i,0]

        d2 = d.reshape(n,n,3)
        plt.contourf(x0[:,:,0], x0[:,:,1], d2[:,:,0])
        cbar = plt.colorbar()
        plt.show()
