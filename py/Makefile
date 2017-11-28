include ../include.mk

alldirs = $(GUI) pyplot

MODULES = __init__.py\
          equi2d.py\
          boundary.py\
          separatrix.py\
          grid.py\
          dataset.py\
          plot_data.py

INSTALLDIR = $(PYTHONDIR)/flare


.PHONY: all install clean uninstall $(alldirs)
all: $(alldirs)


install: $(alldirs)
	install -d $(INSTALLDIR)
	install -m 644 $(MODULES) $(INSTALLDIR)


clean: $(alldirs)


uninstall: $(alldirs)
	-rm -rf $(INSTALLDIR)


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
