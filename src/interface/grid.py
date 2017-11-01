from grid_interface import *

UNSTRUCTURED    = 1
SEMI_STRUCTURED = 2
STRUCTURED      = 3
MESH_2D         = 4
MESH_3D         = 5

CYL_LABEL = ['Major Radius [cm]', 'Vertical coordinate [cm]', 'Toroidal Angle [deg]']

class Grid():
    def __init__(self, filename):
        grid_interface.load(filename)
        self.n = grid_interface.n
        self.x = grid_interface.x
        self.layout = grid_interface.layout

        if self.layout == STRUCTURED:
            self.n1 = grid_interface.n1
            self.n2 = grid_interface.n2
            self.n3 = grid_interface.n3
            self.x1 = grid_interface.x1
            self.x2 = grid_interface.x2
            self.x3 = grid_interface.x3

        self.x1_label = CYL_LABEL[grid_interface.coord1-1]
        self.x2_label = CYL_LABEL[grid_interface.coord2-1]
