include ../../include.mk


OBJECTS =  \
               math.o\
               search.o\
               bilinear_system.o\
               ode_solver.o\
               curve2D.o\
               quad_ele.o\
               usr_grid.o

INC_DIR = -I ..


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
