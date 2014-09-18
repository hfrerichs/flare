#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate structured grid using toroidal angle (phi) and poloidal angle (theta)
#=======================================================================

# select grid position
r_center  = 175.0	# reference position (major radius) for r_min [cm]
r_min     = 47.69	# fixed minor radius [cm]
phi1      = -22.5	# Start position phi [deg]
phi2      = 157.5	#   End position phi [deg]
theta1    = 120.0	# Start position theta [deg]
theta2    = 240.0	#   End position theta [deg]

# select grid resolution
phi_res   = 361		# in phi direction
theta_res = 361		# in theta direction
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 431     (structured grid: poloidal angle theta [deg], toroidal angle phi [deg] at fixed minor radius rmin [cm])\n")
grid_file.printf("# reference radius:  r_center=  %7.2f\n", r_center)
grid_file.printf("# theta resolution:  n_theta =  %8d", theta_res)
grid_file.printf("       theta range:  %6.2f -> %6.2f\n", theta1, theta2)
grid_file.printf("# phi resolution:    n_phi   =  %8d", phi_res)
grid_file.printf("       phi range:    %6.2f -> %6.2f\n", phi1, phi2)
grid_file.printf("# minor radius:      r_min   =  %7.2f\n", r_min)

grid_file.printf("%18.10E\n", theta1)
for j in 2..theta_res
   theta = theta1  +  (theta2-theta1) * (j-1) / (theta_res-1)
   grid_file.printf("%18.10E\n", theta)
end

grid_file.printf("%18.10E\n", phi1)
for i in 2..phi_res
   phi = phi1  +  (phi2-phi1) * (i-1) / (phi_res-1)
   grid_file.printf("%18.10E\n", phi)
end
grid_file.close
#=======================================================================
