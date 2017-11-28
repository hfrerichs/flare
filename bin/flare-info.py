#!/usr/bin/env python
import argparse

from flare.dataset import Dataset



#===============================================================================
# INFO command: query data file for available data
#===============================================================================
def info():
    parser = argparse.ArgumentParser(prog="flare info")
    parser.add_argument("data_file")
    args   = parser.parse_args()
    print Dataset(args.data_file)
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    info()
