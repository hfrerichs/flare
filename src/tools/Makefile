include ../../include.mk


OBJECTS =  \
               trace_bline.o\
               sample_bfield.o\
               poincare_plot.o\
               connection_length.o\
               get_equi_info_2D.o\
               generate_flux_surface_2D.o\
               generate_flux_surface_3D.o\
               plot_boundary.o\
               safety_factor.o\
               transform_to_flux_coordinates.o\
               generate_mag_file.o\
               generate_magnetic_axis.o\
               generate_field_aligned_grid.o\
               flux_surface_grid.o\
               field_line_loss.o\
               generate_separatrix.o\
               footprint_grid.o\
               distance_to_surface.o\
               homoclinic_tangle.o\
               generate_flux_tube.o\
               FLR_analysis.o


INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(BFIELD_DIR) -I ../$(CORE_DIR) -I ../$(GEOMETRY_DIR) -I ../$(GRIDGEN_DIR)


.PHONY: all
all: $(OBJECTS)

%.o: %.f*
	$(FC) -c $(INC_DIR) $<

%.o: %.F*
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
