include ../../include.mk


OBJECTS =  \
               system.o\
               string.o\
               quicksort.o\
               math.o\
               numerics.o\
               run_control.o\
               dataset.o\
               grid.o\
               bspline2D.o\
               cspline.o\
               Q4.o\
               search.o\
               curve2D.o\
               mesh_spacing.o\
               bilinear_system.o\
               interpolate3D.o\
               ode_solver.o\
               quad_ele.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR)


.PHONY: all
all: $(OBJECTS)

%.o: %.f90
	$(FC) -c $(INC_DIR) $<

%.o: %.F90
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
