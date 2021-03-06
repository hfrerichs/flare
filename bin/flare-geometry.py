#!/usr/bin/env PYTHON
import argparse

from flare.backend import grid_interface, mesh_spacing_interface



#===============================================================================
# GRID command: access grid related functions
#===============================================================================

# GRID.1 mesh spacing
def mesh_spacing(args):
    filename = re.sub(" ", "_", args.spacing_command)
    mesh_spacing_interface.generate(args.spacing_command, args.sample_segments, filename, args.output_format)


# GRID.2 generate grid nodes
def nodes(args):
    grid_interface.create(args.n1, args.x1_range[0], args.x1_range[1],
                          args.n2, args.x2_range[0], args.x2_range[1],
                          args.n3, args.x3_range[0], args.x3_range[1],
                          args.output_file)


# GRID.main
def grid():
    # create main parser
    parser     = argparse.ArgumentParser(prog="flare geometry")
    subparsers = parser.add_subparsers(title="commands")


    # GRID.1 mesh/grid spacing
    p_spacing  = subparsers.add_parser("spacing", help="generate mesh/grid spacing function")
    p_spacing.add_argument("spacing_command")
    p_spacing.add_argument("-s", "--sample_segments", default=20, help="number of segments to sample")
    spacing_formats=['function', 'list']
    p_spacing.add_argument("-o", "--output_format", default=spacing_formats[0], choices=spacing_formats)
    p_spacing.set_defaults(func=mesh_spacing)


    # GRID.2 create grid
    p_nodes    = subparsers.add_parser("create", help="create new grid (e.g. for connection length calculation)")
    p_nodes.add_argument("-x1", "--x1_range", nargs=2, type=float, default=[100.0, 160.0])
    p_nodes.add_argument("-n1", default=128, type=int, help="resolution along first coordinate")
    p_nodes.add_argument("-x2", "--x2_range", nargs=2, type=float, default=[-140.0, -80.0])
    p_nodes.add_argument("-n2", default=128, type=int, help="resolution along second coordinate")
    p_nodes.add_argument("-x3", "--x3_range", nargs=2, type=float, default=[0.0, 0.0])
    p_nodes.add_argument("-n3", default=1, type=int, help="resolution along third coordinate")
    p_nodes.add_argument("-o", "--output_file", default="grid.dat")
    p_nodes.set_defaults(func=nodes)


    # parse arguments and execute selected command
    args = parser.parse_args()
    args.func(args)
#===============================================================================




#===============================================================================
if __name__ == "__main__":
    grid()
