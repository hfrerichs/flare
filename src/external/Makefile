include ../../include.mk


OBJECTS =  \
               bspline90_22.o\
               $(ODE_OBJECTS)

EMC3_DIR = EMC3


.PHONY: all
all: $(OBJECTS)
	cd $(EMC3_DIR); $(MAKE) FC='$(FC)'

%.o: %.f*
	$(FC) -c $<

.PHONY: clean
clean:
	cd $(EMC3_DIR); $(MAKE) clean
	rm -rf *.o *.mod
