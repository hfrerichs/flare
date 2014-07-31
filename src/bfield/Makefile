include ../../include.mk


OBJECTS =  \
               reconstruct.o\
               geqdsk.o\
               polygones.o\
               bfield.o\
               equilibrium.o\
               fieldline.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(MATH_DIR) -I ../$(GEOMETRY_DIR)


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
