#!/usr/bin/env python
import argparse
import os, glob

import flare.backend



#===============================================================================
# EQUI3D command: non-axisymmetric (3D) equilibrium related functions
#===============================================================================

# EQUI3D.1 magnetic axis related functions
def equi3d_axis_generate(args):
    x0   = [args.R0, 0.0, 0.0]
    nsym = args.symmetry
    nphi = args.toroidal_resolution
    flare.backend.init(configuration='./')
    flare.backend.generate_magnetic_axis_interface(x0, nsym, nphi)
    for filename in glob.glob('magnetic_axis_*.dat'):
        os.remove(filename)


# EQUI3D.main
def equi3d():
    parser     = argparse.ArgumentParser(prog="flare equi3d")
    subparsers = parser.add_subparsers(title='commands')

    # EQUI3D.1 magnetic axis
    p_axis = subparsers.add_parser("axis", help="magnetic axis related functions")
    s_axis = p_axis.add_subparsers()
    p_axis_generate = s_axis.add_parser("generate", help="generate magnetic axis")
    p_axis_generate.add_argument("R0", type=float, help="Estimated position of magnetic axis at Phi = 0 deg")
    p_axis_generate.add_argument("-s", "--symmetry", default=1, help="toroidal symmetry")
    p_axis_generate.add_argument("-n", "--toroidal-resolution", default=360)
    p_axis_generate.set_defaults(func=equi3d_axis_generate)


    # parse arguments and execute selected command
    args = parser.parse_args()
    args.func(args)
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    equi3d()
