include ../../include.mk

SOURCES = types.f90\
          backend.f90\
          boundary_interface.f90\
          grid_interface.f90\
          equi2d_interface.f90\
          separatrix_interface.f90\
          grid3dgen_interface.f90\
          mesh_spacing_interface.f90

FLARE_SRC_DIR = ../../src
INC_DIRS = "" $(BFIELD_DIR) $(CORE_DIR) $(GEOMETRY_DIR) $(GRIDGEN_DIR)
INC = $(foreach dir,$(INC_DIRS),-I$(FLARE_SRC_DIR)/$(dir))
LIB = -L$(FLARE_SRC_DIR) -lFLARE

BACKEND = backend

INSTALLDIR = $(PYTHONDIR)/flare


.PHONY: all clean install uninstall
all:
	f2py $(INC) -c $(SOURCES) -m $(BACKEND) $(LIB)


clean:
	-rm -f *.so


install:
	install -d $(INSTALLDIR)
	install $(BACKEND).so $(INSTALLDIR)


uninstall:
	-rm -f $(INSTALLDIR)/$(BACKEND).so
