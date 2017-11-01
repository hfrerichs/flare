include ../../include.mk

alldirs = gui

OBJECTS = backend.so\
          grid.so\
          boundary_interface.so\
          python_interface.so

MODULES = __init__.py boundary.py


.PHONY: all clean install uninstall $(alldirs)
all: $(alldirs) $(OBJECTS)


%.so: %.f90
	f2py -I.. -I../$(BFIELD_DIR) -I../$(CORE_DIR) -I../$(GEOMETRY_DIR) -c $< -m $(basename $@) -L.. -lFLARE

clean: $(alldirs)
	-rm -f *.so


install: $(alldirs)
	install -d $(PYTHONDIR)/flare
	install $(OBJECTS) $(PYTHONDIR)/flare
	install -m 644 $(MODULES) $(PYTHONDIR)/flare


uninstall: $(alldirs)
	-rm -rf $(PYTHONDIR)/flare


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
