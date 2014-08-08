include ../include.mk


PROGRAM        = ../bin/flare_bin
PROGRAM_DEBUG  = ../bin/flare_bin_debug

EXTERNAL_INC   = -I $(EXTERNAL_DIR)
EXTERNAL_OBJ   = $(EXTERNAL_DIR)/*.o $(EMC3_OBJ_LONG)

CORE_INC       = -I $(CORE_DIR)
CORE_OBJ       = $(CORE_DIR)/*.o

BFIELD_INC     = -I $(BFIELD_DIR)
BFIELD_OBJ     = $(BFIELD_DIR)/*.o

GEOMETRY_INC   = -I $(GEOMETRY_DIR)
GEOMETRY_OBJ   = $(GEOMETRY_DIR)/*.o

TOOLS_INC      = -I $(TOOLS_DIR)
TOOLS_OBJ      = $(TOOLS_DIR)/*.o

OBJECTS = \
               parallel.o\
               run_control.o


.PHONY: all
all: $(OBJECTS)
	cd $(EXTERNAL_DIR); $(MAKE) FC='$(FC)'
	cd $(CORE_DIR); $(MAKE) FC='$(FC)'
	cd $(BFIELD_DIR); $(MAKE) FC='$(FC)'
	cd $(GEOMETRY_DIR); $(MAKE) FC='$(FC)'
	cd $(TOOLS_DIR); $(MAKE) FC='$(FC)'
	$(FC) -c $(GEOMETRY_INC) $(BFIELD_INC) main.f90
	$(FC) -o  $(PROGRAM) $(OBJECTS) $(TOOLS_OBJ) $(GEOMETRY_OBJ) $(BFIELD_OBJ) $(CORE_OBJ) $(EXTERNAL_OBJ) main.o


%.o: %.f*
	$(FC) -c $<

%.o: %.F*
	$(FC) -c $<


.PHONY: debug
debug:
	$(MAKE) FC="$(FC_DEBUG)" PROGRAM="$(PROGRAM_DEBUG)"

.PHONY: clean
clean:
	cd $(EXTERNAL_DIR); $(MAKE) clean
	cd $(CORE_DIR); $(MAKE) clean
	cd $(GEOMETRY_DIR); $(MAKE) clean
	cd $(BFIELD_DIR); $(MAKE) clean
	cd $(TOOLS_DIR); $(MAKE) clean
	rm -rf *.o
	rm -rf *.mod
