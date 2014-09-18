#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate unstructured grid in RZ-plane (list of R,Z coordinates)
#=======================================================================

# select grid position
r1    = 120.0	# Start position in R [cm]
r2    = 230.0	#   End position in R [cm]
z1    = -55.0	# Start position in Z [cm]
z2    =  55.0	#   End position in Z [cm]
phi0  =   0.0	# Toroidal angle [deg]

# select grid resolution
r_res = 101	# in R direction
z_res = 101	# in Z direction
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
n=r_res*z_res
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 213     (unstructured grid in R,Z plane at phi0 [deg])\n")
grid_file.printf("# resolution:        n_RZ    =  %8d\n", n)
grid_file.printf("# fixed position:    phi0    =  %7.2f\n", phi0)

dr = 0.0
if (r_res > 1)
   dr = (r2-r1) / (r_res-1)
end
dz = 0.0
if (z_res > 1)
   dz = (z2-z1) / (z_res-1)
end

for i in 1..r_res
   for j in 1..z_res
      r =  r1  +  (i-1) * dr
      z =  z1  +  (j-1) * dz
     
      grid_file.printf("%18.10E",   r)
      grid_file.printf("%18.10E\n", z)
   end
end
grid_file.close
#=======================================================================
