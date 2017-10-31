include ../../include.mk

alldirs = gui

MODULES = __init__.py


.PHONY: all clean install uninstall $(alldirs)
all: $(alldirs)
	f2py -I.. -I../$(BFIELD_DIR) -I../$(CORE_DIR) -I../$(GEOMETRY_DIR) -c python_interface.f90 -m backend -L.. -lFLARE


clean: $(alldirs)
	-rm -f *.so


install: $(alldirs)
	install -d $(PYTHONDIR)/flare
	install backend.so $(PYTHONDIR)/flare
	install -m 644 $(MODULES) $(PYTHONDIR)/flare


uninstall: $(alldirs)
	-rm -rf $(PYTHONDIR)/flare


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
