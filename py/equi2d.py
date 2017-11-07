from backend import equi2d_interface


# Load equilibrium data file
# Error codes:
#   1:  invalid format string given
#   2:  no format string given and guess failed
def load(filename, format=None):
    if format:
        iequi = equi2d_interface.get_equilibrium_format_from_string(format)
        if iequi == 0: return 1
    else:
        iequi = equi2d_interface.guess_equilibrium_format(filename)
        if iequi == 0: return 2

    equi2d_interface.load(filename, iequi)

    return 0
