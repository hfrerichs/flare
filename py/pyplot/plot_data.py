import os, sys
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

from flare import Grid, STRUCTURED, SEMI_STRUCTURED


FLARE                = "# FLARE"
FLARE_DATA_DIMENSION = "# FLARE DATA DIMENSION"
FLARE_GEOMETRY       = "# FLARE GEOMETRY"
FLARE_DATA_COLUMN    = "# FLARE DATA COLUMN"
FLARE_DERIVED_DATA   = "# FLARE DERIVED DATA"

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
            if not s.startswith(FLARE): continue

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
                    self.q_derived[d[0]] = (d[1], d[2], d[3])

        f.close()


    def decode_derived_data(self, text):
        # 1. decode data key
        s = text.split()
        if len(s) == 0:
            print "error in data description: no data key given!"
            return None
        q = s[0]


        # 2. decode dependencies
        s = text[len(q):].strip()
        if not s.startswith('['):
            print "error in data description: unexpected character!"
            print "'%s'" % s
            return None

        i = s.find(']')
        if i == -1:
            print "error in data description: unexpected end of dependency list!"
            print "'%s'" % s
            return None
        dependency = [d.strip() for d in s[1:i].split(',')]


        # 3+4. decode recipe and label
        s = s[i+1:].strip()
        s = [d.strip() for d in s.split('"') if d.strip() != '']
        if len(s) < 2:
            print "error in data description: cannot read recipe and label!"
            print "'%s'" %s
            return None

        recipe = s[0]
        for subst in SUBST:
            recipe = re.sub(subst[0], subst[1], recipe)
        label  = s[1]
        return q, dependency, recipe, label


    # returns true if derived data dkey is available
    def supports(self, dkey):
        for key in self.q_derived[dkey][0]:
            if key not in self.q: return False
        return True


    def __str__(self):
        if self.ndim == None  or not self.q:
            return "Data dimension missing or no data description available\n"

        info = "Available data:\n"
        for key, value in self.q.items():
            info += "\t{}: {}\n".format(key, value)

        # Derived data
        n = 0
        derived_info = "Supported derived data:\n"
        for key in self.q_derived:
            if self.supports(key):
                n +=1
                derived_info += "\t{}: {}\n".format(key, self.q_derived[key][2])
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
        for q0 in self.q_derived[qkey][0]:
            exec("{} = self.get_data('{}')".format(q0, q0))
        try:
            exec("d = "+self.q_derived[qkey][1])
        except:
            print "error: cannot evaluate derived quantity {}".format(qkey)
            print "from", self.q_derived[qkey][1]
            sys.exit(2)

        return d


    # return data qkey
    def get_data(self, qkey):
        if qkey in self.q_derived: return self.get_derived_data(qkey)

        i = self.q.keys().index(qkey)
        return self.d[:,i]

    def get_data_label(self, qkey):
        if qkey in self.q_derived: return self.q_derived[qkey][2]

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
    if not qkey in d.available_data():
        print "error: '{}' is not available in data file {}".format(qkey, data_file)
        sys.exit(2)


    # load data
    d.load()



    # select plot type
    ndim = d.get_data_dimension()
    if ndim == PLOT_2D:
        #print "plotting", qkey, "on", grid_file
        d.plot_2d(qkey, grid_file)

    else:
        print "error: plotting of data dimension '{}' is not implemented (yet)!".format(ndim)
        sys.exit(2)
