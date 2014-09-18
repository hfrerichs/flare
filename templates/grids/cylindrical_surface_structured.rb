#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate structured cylindrical surface grid
#=======================================================================

# grid resolution
z_res   = 201		# resolution in vertical direction
phi_res = 241		# resolution in angular/toroidal direction

# grid position
phi1 = -20.0		# Start angle [deg]
phi2 = 100.0		#   End angle [deg]
z1   =-122.3		# Start vertical position [cm]
z2   =-112.3		#   End vertical position [cm]
r    = 101.600000	# Fixed radius [cm]
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 223     (structured cylindrical surface grid)\n")
grid_file.printf("# Z   resolution:    n_Z     =  %8d\n", z_res)
grid_file.printf("# phi resolution:    n_p     =  %8d\n", phi_res)

dz = 0.0
if (z_res > 1)
   dz = (z2-z1) / (z_res-1)
end

for i in 1..z_res
   z = z1  +  dz * (i-1)
   grid_file.printf("%18.10E",   r)
   grid_file.printf("%18.10E",   z)
   grid_file.printf("  %8.3f\n", z)
end

dphi = 0.0
if (phi_res > 0)
   dphi = (phi2-phi1) / (phi_res-1)
end

for j in 1..phi_res
   p = phi1  +  dphi * (j-1)
   grid_file.printf("%18.10E\n", p)
end
grid_file.close
#=======================================================================
