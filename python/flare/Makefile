include ../../include.mk

alldirs = ui

MODULES = __init__.py


.PHONY: all clean install uninstall $(alldirs)
all: $(alldirs)
	f2py -I../../src -I../../src/$(BFIELD_DIR) -I../../src/$(CORE_DIR) -I../../src/$(GEOMETRY_DIR) -c python_interface.f90 -m backend -L../../src -lFLARE


clean: $(alldirs)


install: $(alldirs)
	install -d $(PYTHONDIR)/flare
	install backend.so $(PYTHONDIR)/flare
	install -m 644 $(MODULES) $(PYTHONDIR)/flare


uninstall: $(alldirs)
	-rm -rf $(PYTHONDIR)/flare


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
