import sys

from backend import grid_interface


UNSTRUCTURED    = 1
SEMI_STRUCTURED = 2
STRUCTURED      = 3
MESH_2D         = 4
MESH_3D         = 5
COMPOSITE       = 9

USER_DEFINED    = 0
CARTESIAN       = 1
CYLINDRICAL     = 2

CYL_LABEL = ['Major Radius [cm]', 'Vertical coordinate [cm]', 'Toroidal Angle [deg]']


class Grid():
    def __init__(self, filename):
        if filename.startswith("COMPOSITE"):
            self.load_composite(filename)
        else:
            self.load(filename)


    def load_composite(self, c):
        self.layout = COMPOSITE
        i1 = c.find("(")+1
        i2 = c.rfind(")")
        if i1<0  or  i2<0:
            print("error: invalid geometry '{}'!".format(c))
            sys.exit(1)

        filelist = [s.strip() for s in c[i1:i2].split(',')]
        self.G = []
        for filename in filelist:
            self.G.append(Grid(filename))

        self.n = 0
        for G in self.G:
            self.n += G.n


    def load(self, filename):
        grid_interface.load(filename)
        self.n = grid_interface.n
        self.x = grid_interface.x
        self.layout = grid_interface.layout
        self.coordinates = grid_interface.coordinates

        # set default coordinate labels
        self.x1_label = CYL_LABEL[grid_interface.coord1-1]
        self.x2_label = CYL_LABEL[grid_interface.coord2-1]

        # set dimension of structured grid
        self.n1 = grid_interface.n1
        self.n2 = grid_interface.n2
        self.n3 = grid_interface.n3


        if self.layout == STRUCTURED:
            self.x1 = grid_interface.x1
            self.x2 = grid_interface.x2
            self.x3 = grid_interface.x3


        if self.layout == SEMI_STRUCTURED:
            self.x1 = grid_interface.x1
            self.x2 = grid_interface.x2

            self.x1_label = 'sample coordinate'
            self.x2_label = CYL_LABEL[grid_interface.fixed_coord-1]

        if self.layout == MESH_2D:
            self.x1 = grid_interface.x1
            self.x2 = grid_interface.x2


        if self.coordinates == USER_DEFINED:
            self.x1_label = grid_interface.coord_label1.tostring().strip()
            self.x2_label = grid_interface.coord_label2.tostring().strip()
