include ../../include.mk


OBJECTS =  \
               emc3_grid.o\
               fieldline_grid.o\
               inner_boundary.o\
               divertor.o\
               topo_sc.o\
               topo_lsn.o\
               topo_ddn.o\
               topo_cdn.o\
               base_grid.o\
               mesh_interface.o\
               mfs_mesh.o\
               elements.o\
               base_mesh.o\
               vacuum_domain.o\
               trace_nodes.o\
               check.o


INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(BFIELD_DIR) -I ../$(CORE_DIR) -I ../$(GEOMETRY_DIR)


.PHONY: all
all: $(OBJECTS)

%.o: %.f90
	$(FC) -c $(INC_DIR) $<

%.o: %.F90
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
