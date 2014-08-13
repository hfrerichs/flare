include ../../include.mk


OBJECTS =  \
               boundary.o\
               fieldline.o\
               flux_surface_2D.o

INC_DIR = -I .. -I ../$(CORE_DIR) -I ../$(BFIELD_DIR)


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
