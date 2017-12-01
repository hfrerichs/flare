from copy import copy
import sys
import matplotlib.pyplot as plt
import matplotlib.tri    as tri
from matplotlib.mlab   import griddata as mgriddata
from scipy.interpolate import griddata as sgriddata

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


class Grid(object):
    def __init__(self, filename=None, layout=None):
        self.layout = layout
        if filename:
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
        self.grid = []
        for filename in filelist:
            self.grid.append(Grid(filename))

        self.n = 0
        for grid in self.grid:
            self.n += grid.n


    def load(self, filename):
        grid_interface.load(filename)
        self.n = copy(grid_interface.n)
        self.x = copy(grid_interface.x)
        self.layout = copy(grid_interface.layout)
        self.coordinates = copy(grid_interface.coordinates)

        # set default coordinate labels
        self.x1_label = CYL_LABEL[grid_interface.coord1-1]
        self.x2_label = CYL_LABEL[grid_interface.coord2-1]

        # set dimension of structured grid
        self.n1 = copy(grid_interface.n1)
        self.n2 = copy(grid_interface.n2)
        self.n3 = copy(grid_interface.n3)


        if self.layout == STRUCTURED:
            self.x1 = copy(grid_interface.x1)
            self.x2 = copy(grid_interface.x2)
            self.x3 = copy(grid_interface.x3)


        if self.layout == SEMI_STRUCTURED:
            self.x1 = copy(grid_interface.x1)
            self.x2 = copy(grid_interface.x2)

            self.x1_label = 'sample coordinate'
            self.x2_label = CYL_LABEL[grid_interface.fixed_coord-1]

        if self.layout == MESH_2D:
            self.x1 = copy(grid_interface.x1)
            self.x2 = copy(grid_interface.x2)


        if self.coordinates == USER_DEFINED:
            self.x1_label = grid_interface.coord_label1.tostring().strip()
            self.x2_label = grid_interface.coord_label2.tostring().strip()


    def cells(self):
        if self.layout == COMPOSITE:
            m = 0
            for grid in self.grid:
                m += grid.cells()
            return m

        if self.layout == UNSTRUCTURED:
            return 0

        m = (self.n1-1) * (self.n2-1)
        return m


    def nodes(self):
        return self.n


    # 2d visualization of data in cells
    def plot2d_cells(self, data, *args, **kwargs):
        m = self.cells()
        if len(data) != m:
            print("error: data size ({}) and geometry ({}) are incompatible!".format(len(data), m))
            sys.exit(1)

        if self.layout == COMPOSITE:
            m1 = 0
            for grid in self.grid:
                m  = grid.cells()
                m2 = m1 + m
                print "m1->m2: ", m1, m2
                grid.plot2d_cells(data[m1:m2], *args, **kwargs)
                m1 += m
            return

        if self.layout == MESH_2D:
            x1   = self.x1.reshape(self.n2, self.n1)
            x2   = self.x2.reshape(self.n2, self.n1)
            data = data.reshape(self.n2-1, self.n1-1)
            plt.pcolormesh(x1, x2, data, *args, **kwargs)

        if self.layout == STRUCTURED  or  self.layout == SEMI_STRUCTURED:
            plt.pcolormesh(self.x1, self.x2, data, *args, **kwargs)


    # 2d visualization of data on nodes
    # - unstructured nodes -
    def plot2d_unstructured_nodes(self, data, *args, **kwargs):
        method = kwargs.get('plot_function', 'tricontourf')

        x1   = self.x[:,0]
        x2   = self.x[:,1]
        xmin = min(x1)
        xmax = max(x1)
        ymin = min(x2)
        ymax = max(x2)
        if method == 'tricontourf':
            triang = tri.Triangulation(x1, x2)
            plt.tricontourf(triang, data, *args, **kwargs)

        elif method == 'scipy_griddata_pcontourf':
            n  = np.sqrt(len(data))
            xi = np.linspace(xmin, xmax, n)
            yi = np.linspace(ymin, ymax, n)
            zi = sgriddata((x1, x2), q, (xi[None,:], yi[:,None]), method='nearest')
            cs = plt.contourf(xi, yi, zi, *args, **kwargs)

        elif method == 'matplotlib_griddata_pcontourf':
            n  = np.sqrt(len(data))
            xi = np.linspace(xmin, xmax, n)
            yi = np.linspace(ymin, ymax, n)
            zi = mgriddata(x1, x2, q, xi, yi, interp='linear')
            cs = plt.contourf(xi, yi, zi, *args, **kwargs)


    # - structured nodes -
    def plot2d_structured_nodes(self, data, *args, **kwargs):
        if self.layout == STRUCTURED  or  self.layout == SEMI_STRUCTURED:
            q    = data.reshape(self.n2, self.n1)
            plt.contourf(self.x1, self.x2, q, *args, **kwargs)

        if self.layout == MESH_2D:
            x1   = self.x1.reshape(self.n2, self.n1)
            x2   = self.x2.reshape(self.n2, self.n1)
            q    = q.reshape(self.n2, self.n1)
            plt.contourf(x1, x2, q, *args, **kwargs)


    # - wrapper -
    def plot2d_nodes(self, data, *args, **kwargs):
        if self.layout == UNSTRUCTURED:
            self.plot2d_unstructured_nodes(data, *args, **kwargs)
        else:
            self.plot2d_structured_nodes(data, *args, **kwargs)
