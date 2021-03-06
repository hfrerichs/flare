import re, sys

from dataset import Dataset



# Interface for data visualization
# Input:
#   data_file
#   qkey        data key
#   grid_file   (optional) non-default geometry
def plot_data(data_file, qkey, grid_file=None, dimension=None, *args, **kwargs):
    qkey = re.sub('\$', 'COLUMN', qkey)

    # get data abstract
    d = Dataset(data_file)

    # plot user defined derived quantity
    if re.search('=', qkey):
        recipe = qkey
        qkey   = qkey.split('=')[0]
        d.add_derived_data(qkey, recipe, recipe)

    # plot by column number
    try:
        i = int(qkey)
        recipe = 'X=COLUMN{}'.format(i)
        qkey   = 'X'
        d.add_derived_data(qkey, recipe, "data column {}".format(i))
    except:
        pass


    # check if selected data is available
    if not qkey in d.available_data():
        print "error: '{}' is not available in data file {}".format(qkey, data_file)
        sys.exit(2)


    # load data
    d.load()



    # select plot type
    ndim = d.get_data_dimension()
    if dimension: ndim = dimension
    if ndim == 2:
        #print "plotting", qkey, "on", grid_file
        d.plot_2d(qkey, grid_file, *args, **kwargs)

    else:
        print "error: plotting of data dimension '{}' is not implemented (yet)!".format(ndim)
        sys.exit(2)
