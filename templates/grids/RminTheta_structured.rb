#!/usr/bin/env ruby

require "bigdecimal/math.rb"
pi = Math::PI



#=======================================================================
# Generate structured grid using minor radius (rmin) and poloidal angle (theta)
#=======================================================================

# select grid position
r_center  = 175.0	# reference position (major radius) for rmin [cm]
rmin1     =  38.0	# Start position rmin [cm]
rmin2     =  48.0	#   End position rmin [cm]
theta1    =   0.0	# Start position theta [deg]
theta2    = 360.0	#   End position theta [deg]
phi0      =  0.0	# Toroidal angle [deg]

# select grid resolution
rmin_res  = 101	# in rmin direction
theta_res = 101	# in theta direction
#=======================================================================



#=======================================================================
# Main Program
#=======================================================================
grid_file = File.new("grid.dat", "w")
grid_file.printf("# grid_id = 433     (structured grid: minor radius rmin [cm], poloidal angle theta [deg] at fixed toroidal position phi0 [deg])\n")
grid_file.printf("# reference radius:  r_center=  %7.2f\n", r_center)
grid_file.printf("# rmin resolution:   n_rmin  =  %8d", rmin_res)
grid_file.printf("       rmin range:   %6.2f -> %6.2f\n", rmin1, rmin2)
grid_file.printf("# theta resolution:  n_theta =  %8d", theta_res)
grid_file.printf("       theta range:  %6.2f -> %6.2f\n", theta1, theta2)
grid_file.printf("# fixed position:    phi0    =  %7.2f\n", phi0)

grid_file.printf("%18.10E\n", rmin1)
for i in 2..rmin_res
   rmin = rmin1  +  (rmin2-rmin1) * (i-1) / (rmin_res-1)
   grid_file.printf("%18.10E\n", rmin)
end

grid_file.printf("%18.10E\n", theta1)
for j in 2..theta_res
   theta = theta1  +  (theta2-theta1) * (j-1) / (theta_res-1)
   grid_file.printf("%18.10E\n", theta)
end
grid_file.close
#=======================================================================
