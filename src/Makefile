include ../include.mk


EXTERNAL_INC   = -I $(EXTERNAL_DIR)
EXTERNAL_OBJ   = $(EXTERNAL_DIR)/*.o $(EMC3_OBJ_LONG)

CORE_INC       = -I $(CORE_DIR)
CORE_OBJ       = $(CORE_DIR)/*.o

BFIELD_INC     = -I $(BFIELD_DIR)
BFIELD_OBJ     = $(BFIELD_DIR)/*.o

GEOMETRY_INC   = -I $(GEOMETRY_DIR)
GEOMETRY_OBJ   = $(GEOMETRY_DIR)/*.o

GRIDGEN_INC    = -I $(GRIDGEN_DIR)
GRIDGEN_OBJ    = $(GRIDGEN_DIR)/*.o

TOOLS_INC      = -I $(TOOLS_DIR)
TOOLS_OBJ      = $(TOOLS_DIR)/*.o

DEVEL_INC      = -I $(DEVEL_DIR)
DEVEL_OBJ      = $(DEVEL_DIR)/*.o

OBJECTS = \
               parallel.o


.PHONY: all
all: $(OBJECTS)
	cd $(EXTERNAL_DIR); $(MAKE) FC='$(FC)'
	cd $(CORE_DIR); $(MAKE) FC='$(FC)'
	cd $(BFIELD_DIR); $(MAKE) FC='$(FC)'
	cd $(GEOMETRY_DIR); $(MAKE) FC='$(FC)'
	cd $(GRIDGEN_DIR); $(MAKE) FC='$(FC)'
	cd $(TOOLS_DIR); $(MAKE) FC='$(FC)'
	cd $(DEVEL_DIR); $(MAKE) FC='$(FC)'
	$(FC) -c $(EXTERNAL_INC) $(CORE_INC) $(GEOMETRY_INC) $(BFIELD_INC) main.f90
	$(FC) -o ../bin/$(PROGRAM) $(OBJECTS) $(DEVEL_OBJ) $(TOOLS_OBJ) $(GRIDGEN_OBJ) $(GEOMETRY_OBJ) $(BFIELD_OBJ) $(CORE_OBJ) $(EXTERNAL_OBJ) $(LIBS) main.o


%.o: %.f90
	$(FC) $(CORE_INC) -c $<

%.o: %.F90
	$(FC) $(CORE_INC) -c $<


.PHONY: debug
debug:
	$(MAKE) FC="$(FC_DEBUG)" PROGRAM="$(PROGRAM_DEBUG)"


.PHONY: clean
clean:
	cd $(EXTERNAL_DIR); $(MAKE) clean
	cd $(CORE_DIR); $(MAKE) clean
	cd $(GEOMETRY_DIR); $(MAKE) clean
	cd $(BFIELD_DIR); $(MAKE) clean
	cd $(GRIDGEN_DIR); $(MAKE) clean
	cd $(TOOLS_DIR); $(MAKE) clean
	cd $(DEVEL_DIR); $(MAKE) clean
	rm -rf *.o
	rm -rf *.mod


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd $(ADDONS_DIR); $(MAKE) $(MAKECMDGOALS)
