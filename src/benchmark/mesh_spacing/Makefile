include ../../../include.mk

PROGRAM = benchmark_mesh_spacings

CORE_OBJECTS = search.o math.o quicksort.o dataset.o curve2D.o string.o mesh_spacing.o
CORE_OBJECTS_LONG = $(addprefix ../../$(CORE_DIR)/, $(CORE_OBJECTS))
FLARE_OBJECTS = parallel.o
FLARE_OBJECTS_LONG = $(addprefix ../../, $(FLARE_OBJECTS))
SOURCE = benchmark_mesh_spacings.f90


INC_DIR = -I ../../$(CORE_DIR)


.PHONY: all
all:
	$(FC) $(INC_DIR) $(SOURCE) -o $(PROGRAM) $(FLARE_OBJECTS_LONG) $(CORE_OBJECTS_LONG)


.PHONY: clean
clean:
	rm -rf *.o *.mod *.plt
