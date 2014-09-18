#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate line segment xi->xd with n sample points
#=======================================================================

# initial location (e.g. camera position)
xi = -172.45
yi =  222.32
zi =    0.0

# destination (e.g. observation point)
xd =  54.403
yd = 115.09
zd =   0.0

# relative coordinate of final position
scale = 2.0

# grid resolution
n  = 801
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 110     (Line segment, Cartesian coordinates)\n")
grid_file.printf("# total resolution:  n_xyz   =  %8d\n", n)
for i in 1..n
    x = xi  +  (xd-xi) * (i-1) / (n-1) * scale
    y = yi  +  (yd-yi) * (i-1) / (n-1) * scale
    z = zi  +  (zd-zi) * (i-1) / (n-1) * scale
     
    grid_file.printf("%18.10E",   x)
    grid_file.printf("%18.10E",   y)
    grid_file.printf("%18.10E\n", z)
end
grid_file.close
