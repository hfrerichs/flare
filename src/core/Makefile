include ../../include.mk


OBJECTS =  \
               math.o\
               ode_solver.o


.PHONY: all
all: $(OBJECTS)

%.o: %.f*
	$(FC) -c $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
