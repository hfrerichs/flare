include ../../include.mk


OBJECTS =  \
               reconstruct.o\
               geqdsk.o\
               equilibrium.o\
               polygones.o\
               bfield.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(CORE_DIR)


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
