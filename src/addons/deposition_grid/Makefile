include ../../../include.mk


PROGRAM = deposition_grid

OBJECTS = deposition_grid.o

FLARE_OBJ = parallel.o\
            math.o\
            grid.o


INC_DIR = -I ../.. -I ../../$(CORE_DIR) -I ../../$(GEOMETRY_DIR) -I ../../$(GRIDGEN_DIR)


.PHONY: all clean install
all: $(FLARE_OBJ) $(OBJECTS)
	$(FC_ADDONS) -o $(PROGRAM) $(OBJECTS) $(FLARE_OBJ)


%.o: %.f90
	$(FC_ADDONS) -cpp -c $(INC_DIR) $<


parallel.f90:
	ln -sf ../../parallel.F90 parallel.f90
math.f90:
	ln -sf ../../$(CORE_DIR)/math.f90
grid.f90:
	ln -sf ../../$(CORE_DIR)/grid.F90 grid.f90


install:
	cp $(PROGRAM) $(BINDIR)


clean:
	rm -rf *.o
	rm -rf *.mod
	rm -f parallel.* math.* grid.*
