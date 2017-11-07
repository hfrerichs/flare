import sys

from backend import equi2d_interface


LOAD_ERROR = ['',
    'invalid format string given',
    'no format string given and guess failed',
    'loading data file failed']


class LoadError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def init_config(filename, format):
    ierr = equi2d_interface.init_config(filename, format)
    if ierr > 0:
        raise LoadError(LOAD_ERROR[ierr])
