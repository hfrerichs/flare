include ../../include.mk

alldirs = gui

SOURCES = types.f90\
          backend.f90\
          boundary_interface.f90\
          grid3dgen_interface.f90

BACKEND = backend
MODULES = __init__.py grid.py boundary.py


.PHONY: all clean install uninstall $(alldirs)
all: $(alldirs)
	f2py -I.. -I../$(BFIELD_DIR) -I../$(CORE_DIR) -I../$(GEOMETRY_DIR) -I../$(GRIDGEN_DIR) -c $(SOURCES) -m $(BACKEND) -L.. -lFLARE


clean: $(alldirs)
	-rm -f *.so


install: $(alldirs)
	install -d $(PYTHONDIR)/flare
	install $(BACKEND).so $(PYTHONDIR)/flare
	install -m 644 $(MODULES) $(PYTHONDIR)/flare


uninstall: $(alldirs)
	-rm -rf $(PYTHONDIR)/flare


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
