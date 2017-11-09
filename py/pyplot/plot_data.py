import os, sys
from collections import OrderedDict


FLARE             = "# FLARE"
FLARE_DATA_TYPE   = "# FLARE DATA TYPE"
FLARE_GEOMETRY    = "# FLARE GEOMETRY"
FLARE_DATA_COLUMN = "# FLARE DATA COLUMN"


class Abstract():
    def __init__(self, data_file):
        self.q         = OrderedDict()
        self.data_type = None

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


    def available_data(self):
        return self.q.keys()


    def get_column_number(self, key):
        return self.q.keys().index(key)



# Interface for data visualization
# Input:
#   data_file
#   q           data key
#   grid_file   (optional) non-default geometry
def plot_data(data_file, q, grid_file=None):
    pass
