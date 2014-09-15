include ../../include.mk


OBJECTS =  \
               block_limiter.o\
               boundary.o\
               fieldline.o\
               poincare_set.o\
               flux_surface_2D.o\
               flux_surface_3D.o\
               separatrix.o\
               field_aligned_grid.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(CORE_DIR) -I ../$(BFIELD_DIR)


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
