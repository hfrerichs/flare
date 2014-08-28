include ../../include.mk


OBJECTS =  \
               hyperbolic_fixed_point.o\
               run_control_development.o


INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(BFIELD_DIR) -I ../$(CORE_DIR) -I ../$(GEOMETRY_DIR)


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
