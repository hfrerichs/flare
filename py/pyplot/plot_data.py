import os, sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri    as tri
from matplotlib.mlab   import griddata as mgriddata
from scipy.interpolate import griddata as sgriddata
from collections       import OrderedDict

from flare import Grid, UNSTRUCTURED, STRUCTURED, SEMI_STRUCTURED, MESH_2D


FLARE_DATA_DIMENSION = "# DATA DIMENSION"
FLARE_GEOMETRY       = "# GEOMETRY"
FLARE_DATA_COLUMN    = "# DATA COLUMN"
FLARE_DERIVED_DATA   = "# DERIVED DATA"

PLOT_2D = "2D"

SUBST = [
    ['min', 'np.minimum'],
    ['max', 'np.maximum']
]


class Data():
    def __init__(self, data_file):
        self.q         = OrderedDict()
        self.q_derived = OrderedDict()
        self.ndim      = None
        self.data_file = data_file
        self.geometry  = None

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
                self.geometry = s[len(FLARE_GEOMETRY):].strip()
                #print "geometry = '{}'".format(self.geometry)

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


    # returns true if derived data dkey is available
    def __str__(self):
        if self.ndim == None  or not self.q:
            return "Data dimension missing or no data description available\n"

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
        return self.ndim


    def get_geometry(self):
        return self.geometry


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
        for q0 in self.q:
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
        if qkey in self.q_derived: return self.get_derived_data(qkey)

        i = self.q.keys().index(qkey)
        return self.d[:,i]

    def get_data_label(self, qkey):
        if qkey in self.q_derived: return self.q_derived[qkey][1]

        return self.q[qkey][1:-1]


    def plot_2d(self, qkey, geometry=None, *args, **kwargs):
        # load user defined geometry
        if geometry:
            G = Grid(geometry)
        # load default geometry
        else:
            if not self.geometry:
                print "error: geometry is undefined!"
                sys.exit(2)
            G = Grid(self.geometry)


        # prepare data
        q    = self.get_data(qkey)
        qmin = np.nanmin(q)
        qmax = np.nanmax(q)


        # create plot
        levels = np.linspace(qmin, qmax, 64)
        if G.layout == STRUCTURED  or  G.layout == SEMI_STRUCTURED:
            q    = q.reshape(G.n2, G.n1)
            plt.contourf(G.x1, G.x2, q, vmin=qmin, vmax=qmax, levels=levels, *args, **kwargs)

        if G.layout == MESH_2D:
            x1   = G.x1.reshape(G.n2, G.n1)
            x2   = G.x2.reshape(G.n2, G.n1)
            q    = q.reshape(G.n2, G.n1)
            plt.contourf(x1, x2, q, vmin=qmin, vmax=qmax, levels=levels, *args, **kwargs)

        elif G.layout == UNSTRUCTURED:
            method = kwargs.get('plot_function', 'tricontourf')

            x1   = G.x[:,0]
            x2   = G.x[:,1]
            xmin = min(x1)
            xmax = max(x1)
            ymin = min(x2)
            ymax = max(x2)
            if method == 'tricontourf':
                triang = tri.Triangulation(x1, x2)
                plt.tricontourf(triang, q, levels=levels, *args, **kwargs)

            elif method == 'scipy_griddata_pcontourf':
                n  = np.sqrt(len(q))
                xi = np.linspace(xmin, xmax, n)
                yi = np.linspace(ymin, ymax, n)
                zi = sgriddata((x1, x2), q, (xi[None,:], yi[:,None]), method='nearest')
                cs = plt.contourf(xi, yi, zi, levels=levels, *args, **kwargs)

            elif method == 'matplotlib_griddata_pcontourf':
                n  = np.sqrt(len(q))
                xi = np.linspace(xmin, xmax, n)
                yi = np.linspace(ymin, ymax, n)
                zi = mgriddata(x1, x2, q, xi, yi, interp='linear')
                cs = plt.contourf(xi, yi, zi, levels=levels, *args, **kwargs)

        plt.xlabel(G.x1_label)
        plt.ylabel(G.x2_label)
        cbar = plt.colorbar()
        qlabel = self.get_data_label(qkey)
        cbar.set_label(qlabel)



# Interface for data visualization
# Input:
#   data_file
#   qkey        data key
#   grid_file   (optional) non-default geometry
def plot_data(data_file, qkey, grid_file=None, *args, **kwargs):
    # get data abstract
    d = Data(data_file)
    if not qkey in d.available_data():
        print "error: '{}' is not available in data file {}".format(qkey, data_file)
        sys.exit(2)


    # load data
    d.load()



    # select plot type
    ndim = d.get_data_dimension()
    if ndim == PLOT_2D:
        #print "plotting", qkey, "on", grid_file
        d.plot_2d(qkey, grid_file, *args, **kwargs)

    else:
        print "error: plotting of data dimension '{}' is not implemented (yet)!".format(ndim)
        sys.exit(2)
