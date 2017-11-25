include ../../include.mk


OBJECTS =  \
               run_control_development.o\
               debug.o


INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(BFIELD_DIR) -I ../$(CORE_DIR) -I ../$(GEOMETRY_DIR) -I ../$(GRIDGEN_DIR)


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
