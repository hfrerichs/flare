include ../../../include.mk

PROGRAM = benchmark_mesh_spacing

CORE_OBJECTS = search.o math.o quicksort.o dataset.o curve2D.o string.o mesh_spacing.o
CORE_OBJECTS_LONG = $(addprefix ../../$(CORE_DIR)/, $(CORE_OBJECTS))
FLARE_OBJECTS = parallel.o
FLARE_OBJECTS_LONG = $(addprefix ../../, $(FLARE_OBJECTS))
OBJECTS = benchmark_mesh_spacing.o run_benchmark.o


INC_DIR = -I ../../$(CORE_DIR)


.PHONY: all
all: $(OBJECTS)
	$(FC) -o $(PROGRAM) $(OBJECTS) $(FLARE_OBJECTS_LONG) $(CORE_OBJECTS_LONG)


%.o: %.f90
	$(FC) -c $(INC_DIR) $<


.PHONY: clean
clean:
	rm -rf *.o *.mod *.plt
