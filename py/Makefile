include ../include.mk

alldirs = gui pyplot

MODULES = __init__.py\
          boundary.py\
          grid.py

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
