include ../../include.mk


OBJECTS =  \
               system.o\
               exceptions.o\
               string.o\
               quicksort.o\
               math.o\
               numerics.o\
               linalg.o\
               diffusion.o\
               run_control.o\
               dataset.o\
               grid.o\
               bspline1D.o\
               bspline2D.o\
               bspline3D.o\
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
