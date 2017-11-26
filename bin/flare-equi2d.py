#!/usr/bin/env python
import argparse

import flare.equi2d



#===============================================================================
# EQUI2D command: access axisymmetric (2D) equilibrium data
#===============================================================================

# EQUI2D.1 visualization
def equi2d_view(args):
    flare.equi2d.view(args.data_file, args.format)


# EQUI2D.2 create FLARE configuration for new equilibrium
def equi2d_init(args):
    flare.equi2d.init(args.data_file, args.format)


# EQUI2D.main
def equi2d():
    parser     = argparse.ArgumentParser(prog="flare equi2d")
    subparsers = parser.add_subparsers(title='commands')

    # EQUI2D.1 visualize equilibrium
    p_view = subparsers.add_parser("view", help="visualize equilibrium")
    p_view.add_argument("data_file")
    p_view.add_argument("-f", "--format", default="", help="equilibrium format string")
    p_view.set_defaults(func=equi2d_view)


    # EQUI2D.2 create FLARE configuration
    p_view = subparsers.add_parser("init", help="initialize FLARE configuration for equilibrium")
    p_view.add_argument("-d", "--data_file", help="equilibrium data file (default: use bfield.conf)")
    p_view.add_argument("-f", "--format", default="", help="equilibrium format string")
    p_view.set_defaults(func=equi2d_init)


    # parse arguments and execute selected command
    args = parser.parse_args()
    args.func(args)
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    equi2d()
