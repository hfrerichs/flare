import os, sys
import re
import numpy as np
import matplotlib.pyplot as plt
from collections       import OrderedDict

from flare import Grid, UNSTRUCTURED, STRUCTURED, SEMI_STRUCTURED, MESH_2D


FLARE_DATA_DIMENSION = "# DATA DIMENSION"
FLARE_GEOMETRY       = "# GEOMETRY"
FLARE_DATA_TYPE      = "# TYPE"
FLARE_DATA_COLUMN    = "# DATA COLUMN"
FLARE_DERIVED_DATA   = "# DERIVED DATA"

PLOT_2D = "2D"
FILE       = "FILE"
POINT_DATA = "POINT_DATA"
CELL_DATA  = "CELL_DATA"
NODE_DATA  = "NODE_DATA"

SUBST = [
    ['min', 'np.minimum'],
    ['max', 'np.maximum']
]


class Dataset():
    def __init__(self, data_file):
        self.q         = OrderedDict()
        self.q_derived = OrderedDict()
        self.ndim      = None
        self.data_file = data_file
        self.data_type = NODE_DATA
        self.geometry  = OrderedDict()

        if not os.path.isfile(data_file):
            print "error: data file '{}' does not exist!".format(data_file)
            sys.exit(2)

        f = open(data_file)
        while True:
            s = f.readline()
            if not s: break
            if not s.startswith('# '): break

            # set data dimension
            if s.startswith(FLARE_DATA_DIMENSION):
                self.ndim = s[len(FLARE_DATA_DIMENSION):].strip()
                #print "data dimension = '{}'".format(self.ndim)

            # set geometry
            if s.startswith(FLARE_GEOMETRY):
                g = s[len(FLARE_GEOMETRY):].strip()
                if g.startswith(FILE):
                    filename = g[len(FILE):].strip()
                    self.geometry[filename] = FILE
                else:
                    self.geometry[g] = g

            # set data type
            if s.startswith(FLARE_DATA_TYPE):
                self.data_type = s[len(FLARE_DATA_TYPE):].strip()

            # add data column
            if s.startswith(FLARE_DATA_COLUMN):
                c = s[len(FLARE_DATA_COLUMN):].strip()
                if c == "":
                    print "error: missing data description!"
                    print s
                    sys.exit(2)

                data_id     = c.split()[0]
                description = c[len(data_id):].strip()
                self.q[data_id] = description
                #print "'{}': '{}'".format(data_id, description)

            # add derived data
            if s.startswith(FLARE_DERIVED_DATA):
                c = s[len(FLARE_DERIVED_DATA):].strip()
                d = self.decode_derived_data(c)
                if d:
                    self.q_derived[d[0]] = (d[1], d[2])

        f.close()


    def decode_derived_data(self, text):
        # 1. decode data key
        s = text.split('=')
        if len(s) == 0:
            print "error in data description: no data key given!"
            return None
        q = s[0]


        # 2+3. decode recipe and label
        s = text[len(q)+1:]
        s = [d.strip() for d in s.split('"') if d.strip() != '']
        if len(s) < 2:
            print "error in data description: cannot read recipe and label!"
            print "'%s'" %s
            return None

        recipe = s[0]
        for subst in SUBST:
            recipe = re.sub(subst[0], subst[1], recipe)
        label  = s[1]
        return q, recipe, label


    def add_derived_data(self, dkey, recipe, label):
        if dkey in self.q_derived:
            print "error: {} is already defined!"
            sys.exit(2)

        for subst in SUBST:
            recipe = re.sub(subst[0], subst[1], recipe)
        self.q_derived[dkey] = (recipe, label)


    # returns true if derived data dkey is available
    def __str__(self):
        if not self.q:
            return "no data description available\n"

        info = "Available data:\n"
        for key, value in self.q.items():
            info += "\t{}: {}\n".format(key, value)

        # Derived data
        n = 0
        derived_info = "Supports derived data:\n"
        for key in self.q_derived:
            n +=1
            derived_info += "\t{}: {}\n".format(key, self.q_derived[key][1])
        if n > 0: info += derived_info

        return info


    def get_data_dimension(self):
        ndim = {
            PLOT_2D: 2,
        }.get(self.ndim)
        return ndim


    def make_geometry(self, recipe):
        xy = recipe.split('-')
        geometry = Grid(layout=UNSTRUCTURED)
        x1       = self.get_data(xy[0])
        x2       = self.get_data(xy[1])
        geometry.x = np.concatenate((x1, x2)).reshape(2, len(x1)).T
        geometry.x1_label = self.get_data_label(xy[0])
        geometry.x2_label = self.get_data_label(xy[1])
        return geometry


    def get_geometry(self, user_geometry=None):
        if user_geometry:
            if user_geometry in self.geometry:
                recipe = self.geometry[user_geometry]
                return self.make_geometry(recipe)
            else:
                return Grid(user_geometry)

        # try default geometry
        if len(self.geometry) == 0:
            print "error: geometry is undefined!"
            sys.exit(1)
        key    = self.geometry.keys()[0]
        recipe = self.geometry[key]
        if recipe == FILE:
            return Grid(key)
        else:
            return self.make_geometry(recipe)


    def available_data(self):
        return self.q.keys() + self.q_derived.keys()


    def get_column_number(self, key):
        return self.q.keys().index(key)


    # load data
    def load(self):
        self.d = np.loadtxt(self.data_file, dtype='float')
        if self.d.ndim == 1: self.d = self.d.reshape((self.d.size, 1))
        if self.d.ndim != 2:
            print "error: unexpected dimension of data array!"
            print "dimension = ", self.d.ndim
            sys.exit(2)


    # return derived data
    def get_derived_data(self, qkey):
        # prepare labeled data (by key)
        for q0 in self.q:
            if re.search(q0, self.q_derived[qkey][0]):
                exec("{} = self.get_data('{}')".format(q0, q0))
        # prepare raw data (by column number)
        for i in range(self.d.shape[1]):
            q0 = 'COLUMN{}'.format(i+1)
            if re.search(q0, self.q_derived[qkey][0]):
                exec("{} = self.get_data('{}')".format(q0, q0))

        try:
            exec("d = "+self.q_derived[qkey][0])
        except:
            print "error: cannot evaluate derived quantity {}".format(qkey)
            print "from", self.q_derived[qkey][0]
            sys.exit(2)

        return d


    # return data qkey
    def get_data(self, qkey):
        # 1. qkey is derived data
        if qkey in self.q_derived: return self.get_derived_data(qkey)
        # 2. qkey is primary data
        if qkey in self.q:
            i = self.q.keys().index(qkey)
            return self.d[:,i]
        # 3. qkey is column number
        if qkey.startswith('COLUMN'):
            i = int(qkey[6:])
            return self.d[:,i-1]

    def get_data_label(self, qkey):
        if qkey in self.q_derived: return self.q_derived[qkey][1]

        return self.q[qkey][1:-1]


    def plot_2d(self, qkey, geometry=None, *args, **kwargs):
        geometry = self.get_geometry(user_geometry=geometry)

        # prepare data
        q    = self.get_data(qkey)
        qmin = np.nanmin(q)
        qmax = np.nanmax(q)


        # create plot
        levels = np.linspace(qmin, qmax, 64)

        if self.data_type == CELL_DATA:
            geometry.plot2d_cells(q, vmin=qmin, vmax=qmax)
        elif self.data_type == NODE_DATA:
            geometry.plot2d_nodes(q, vmin=qmin, vmax=qmax, levels=levels)
        elif self.data_type == POINT_DATA:
            geometry.plot2d_points(c=q)
        else:
            print("error: invalid data type '{}'!".format(data_type))
            sys.exit(1)


        plt.xlabel(geometry.x1_label)
        plt.ylabel(geometry.x2_label)
        cbar = plt.colorbar()
        qlabel = self.get_data_label(qkey)
        cbar.set_label(qlabel)
