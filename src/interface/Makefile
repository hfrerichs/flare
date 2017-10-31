include ../../include.mk

MODULES = __init__.py


.PHONY: all clean install uninstall $(alldirs)
all:
	f2py -I.. -I../$(BFIELD_DIR) -I../$(CORE_DIR) -I../$(GEOMETRY_DIR) -c python_interface.f90 -m backend -L.. -lFLARE


clean:
	-rm -f *.so


install:
	install -d $(PYTHONDIR)/flare
	install backend.so $(PYTHONDIR)/flare
	install -m 644 $(MODULES) $(PYTHONDIR)/flare


uninstall:
	-rm -rf $(PYTHONDIR)/flare
