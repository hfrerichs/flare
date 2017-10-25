from PySide.QtGui import *
from PySide.QtCore import *

import flare
from main_ui import Ui_MainWindow


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
        flare.backend.trace_bline_py(a, 1000)
