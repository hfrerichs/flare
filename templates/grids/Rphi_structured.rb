#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate structured cylindrical surface grid
#=======================================================================

# grid resolution
r_res   = 201		# resolution in horizontal direction
phi_res = 241		# resolution in angular/toroidal direction

# grid position
phi1 = -20.0		# Start angle [deg]
phi2 = 100.0		#   End angle [deg]
r1   = 120.0		# Start horizontal position [cm]
r2   = 140.0		#   End horizontal position [cm]
z0   = -123.000		# Fixed vertical position [cm]
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 223     (structured cylindrical surface grid)\n")
grid_file.printf("# R   resolution:    n_Z     =  %8d\n", r_res)
grid_file.printf("# phi resolution:    n_p     =  %8d\n", phi_res)

dr = 0.0
if (r_res > 1)
   dr = (r2-r1) / (r_res-1)
end

for i in 1..r_res
   r = r1  +  dr * (i-1)
   grid_file.printf("%18.10E",   r)
   grid_file.printf("%18.10E",   z0)
   grid_file.printf("  %8.3f\n", r)
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
