from backend import grid_interface

UNSTRUCTURED    = 1
SEMI_STRUCTURED = 2
STRUCTURED      = 3
MESH_2D         = 4
MESH_3D         = 5

USER_DEFINED    = 0
CARTESIAN       = 1
CYLINDRICAL     = 2

CYL_LABEL = ['Major Radius [cm]', 'Vertical coordinate [cm]', 'Toroidal Angle [deg]']


class Grid():
    def __init__(self, filename):
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
