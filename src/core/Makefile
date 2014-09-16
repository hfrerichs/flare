include ../../include.mk


OBJECTS =  \
               math.o\
               dataset.o\
               search.o\
               mesh_spacing.o\
               bilinear_system.o\
               interpolate3D.o\
               ode_solver.o\
               curve2D.o\
               quad_ele.o\
               usr_grid.o\
               grid.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR)


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
