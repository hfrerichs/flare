include ../../include.mk


OBJECTS =  \
               block_limiter.o\
               boundary.o\
               fieldline.o\
               poincare_set.o\
               flux_surface_2D.o\
               flux_surface_3D.o\
               separatrix.o\
               xpaths.o\
               flux_tube.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(CORE_DIR) -I ../$(BFIELD_DIR)


.PHONY: all
all: $(OBJECTS)

%.o: %.f
	$(FC) -c $(INC_DIR) $<

%.o: %.f90
	$(FC) -c $(INC_DIR) $<

%.o: %.F90
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
