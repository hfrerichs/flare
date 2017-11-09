import os, sys
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

from flare import Grid, STRUCTURED, SEMI_STRUCTURED


FLARE             = "# FLARE"
FLARE_DATA_TYPE   = "# FLARE DATA TYPE"
FLARE_GEOMETRY    = "# FLARE GEOMETRY"
FLARE_DATA_COLUMN = "# FLARE DATA COLUMN"

PLOT_2D = "2D"

DERIVED_DATA = ['Lc', 'Lcs']
#    'Lc':   'abs(Lc_bwd) + abs(Lc_fwd)'
#    'Lcs':  'abs(Lc_bwd) + abs(Lc_fwd)'
#}

class Data():
    def __init__(self, data_file):
        self.q         = OrderedDict()
        self.data_type = None
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
            if not s.startswith(FLARE): continue

            # set data type
            if s.startswith(FLARE_DATA_TYPE):
                self.data_type = s[len(FLARE_DATA_TYPE):].strip()
                #print "data type = '{}'".format(self.data_type)

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

        f.close()


    def __str__(self):
        if self.data_type == None  or not self.q:
            return "Data type missing or no data description available\n"

        info = "Available data:\n"
        for key, value in self.q.items():
            info += "\t{}: {}\n".format(key, value)
        return info


    def get_data_type(self):
        return self.data_type


    def get_geometry(self):
        return self.geometry


    def available_data(self):
        return self.q.keys()


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
        d = {
            'Lc':   abs(self.get_data('Lc_bwd'))+abs(self.get_data('Lc_fwd')),
            'Lcs':  np.minimum(abs(self.get_data('Lc_bwd')), abs(self.get_data('Lc_fwd')))
        }.get(qkey)
        return d


    # return data qkey
    def get_data(self, qkey):
        if qkey in DERIVED_DATA: return self.get_derived_data(qkey)

        i = self.q.keys().index(qkey)
        return self.d[:,i]

    def get_data_label(self, qkey):
        if qkey in DERIVED_DATA: return "derived data"

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
        # check geometry
        if G.layout != STRUCTURED  and  G.layout != SEMI_STRUCTURED:
            print "error: visualization on unstructured grids not implemented yet!"
            sys.exit(2)


        # prepare data
        q    = self.get_data(qkey).reshape(G.n2, G.n1)
        qmin = np.nanmin(q)
        qmax = np.nanmax(q)


        # create plot
        levels = np.linspace(qmin, qmax, 64)
        plt.contourf(G.x1, G.x2, q, vmin=qmin, vmax=qmax, levels=levels, *args, **kwargs)
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
def plot_data(data_file, qkey, grid_file=None):
    # get data abstract
    d = Data(data_file)
    if not qkey in d.available_data()  and  not qkey in DERIVED_DATA:
        print "error: '{}' is not available in data file {}".format(qkey, data_file)
        sys.exit(2)


    # load data
    d.load()



    # select plot type
    data_type = d.get_data_type()
    if data_type == PLOT_2D:
        #print "plotting", qkey, "on", grid_file
        d.plot_2d(qkey, grid_file)

    else:
        print "error: plotting of data type '{}' is not implemented (yet)!".format(data_type)
        sys.exit(2)
