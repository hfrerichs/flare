include ../../include.mk


OBJECTS =  \
               pfc.o

INC_DIR = -I .. -I ../$(GEOMETRY_DIR) -I ../$(MATH_DIR) -I ../$(BFIELD_DIR)


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
