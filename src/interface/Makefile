include ../../include.mk

SOURCES = types.f90\
          backend.f90\
          boundary_interface.f90\
          grid_interface.f90\
          equi2d_interface.f90\
          separatrix_interface.f90\
          grid3dgen_interface.f90\
          mesh_spacing_interface.f90

BACKEND = backend

INSTALLDIR = $(PYTHONDIR)/flare


.PHONY: all clean install uninstall
all:
	f2py -I.. -I../$(BFIELD_DIR) -I../$(CORE_DIR) -I../$(GEOMETRY_DIR) -I../$(GRIDGEN_DIR) -c $(SOURCES) -m $(BACKEND) -L.. -lFLARE


clean:
	-rm -f *.so


install:
	install -d $(INSTALLDIR)
	install $(BACKEND).so $(INSTALLDIR)


uninstall:
	-rm -f $(INSTALLDIR)/$(BACKEND).so
