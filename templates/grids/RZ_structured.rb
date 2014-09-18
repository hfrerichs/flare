#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate structured grid in RZ-plane (separate list of R coordinates
# and list of Z coordinates)
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
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 233     (structured grid in R,Z plane at phi0 [deg])\n")
grid_file.printf("# R resolution:      n_R     =  %8d", r_res)
grid_file.printf("           R range:  %6.2f -> %6.2f\n", r1, r2)
grid_file.printf("# Z resolution:      n_Z     =  %8d", z_res)
grid_file.printf("           Z range:  %6.2f -> %6.2f\n", z1, z2)
grid_file.printf("# Fixed position:    phi0    =  %7.2f\n", phi0)

dr = 0.0
if (r_res > 1)
   dr = (r2-r1) / (r_res-1)
end
dz = 0.0
if (z_res > 1)
   dz = (z2-z1) / (z_res-1)
end


for i in 1..r_res
   r = r1  +  (i-1) * dr
   grid_file.printf("%18.10E\n", r)
end

for j in 1..z_res
   z = z1  +  (j-1) * dz
   grid_file.printf("%18.10E\n", z)
end
grid_file.close
#=======================================================================
