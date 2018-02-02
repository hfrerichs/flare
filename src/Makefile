include ../include.mk

alldirs = $(EXTERNAL_DIR) $(CORE_DIR) $(BFIELD_DIR) $(GEOMETRY_DIR) $(GRIDGEN_DIR) $(TOOLS_DIR) $(DEVEL_DIR)

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

VERSION_F90    = version.f90

OBJECTS = \
               version.o\
               parallel.o


.PHONY: all
all: $(OBJECTS) $(alldirs) $(FLARELIB)
	$(FC) -c $(EXTERNAL_INC) $(CORE_INC) $(GEOMETRY_INC) $(BFIELD_INC) main.f90
	$(FC) -o ../bin/$(PROGRAM) main.o $(FLARELIB) $(FIO_PATH)


.PHONY: $(VERSION_F90)
$(VERSION_F90):
	rm -f $(VERSION_F90)
	echo "subroutine print_version()" > $(VERSION_F90)
	echo "  write(6, *) 'Running FLARE'" >> $(VERSION_F90)
	echo -n "  write(6, *) ' " >> $(VERSION_F90)
	git log --pretty=format:"%d" -1 | cut -c 1-72 | tr -d "\n" >> $(VERSION_F90)
	echo "'" >> $(VERSION_F90)
	echo -n "  write(6, *) '  " >> $(VERSION_F90)
	git log --pretty=format:"%ci %h" -1 >> $(VERSION_F90)
	echo "'" >> $(VERSION_F90)
	echo "end subroutine print_version" >> $(VERSION_F90)


.PHONY: $(FLARELIB)
$(FLARELIB):
	$(FC) -shared -o $(FLARELIB) $(OBJECTS) $(DEVEL_OBJ) $(TOOLS_OBJ) $(GRIDGEN_OBJ) $(GEOMETRY_OBJ) $(BFIELD_OBJ) $(CORE_OBJ) $(EXTERNAL_OBJ) $(LIBS)


%.o: %.f90
	$(FC) $(CORE_INC) -c $<

%.o: %.F90
	$(FC) $(CORE_INC) -c $<


.PHONY: debug
debug:
	$(MAKE) FC="$(FC_DEBUG)"


.PHONY: install
install:
	install -d $(LIBDIR)
	install $(FLARELIB) $(LIBDIR)


.PHONY: uninstall
uninstall:
	-rm -f $(LIBDIR)/$(FLARELIB)
	-rm -d $(LIBDIR)


.PHONY: clean
clean: $(alldirs)
	rm -rf *.o
	rm -rf *.mod


.PHONY: $(alldirs)
$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd $(ADDONS_DIR); $(MAKE) $(MAKECMDGOALS)
