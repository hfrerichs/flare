include ../include.mk

MODULES = __init__.py\
          boundary.py\
          grid.py

INSTALLDIR = $(PYTHONDIR)/flare


.PHONY: all install clean uninstall
all:


install:
	install -d $(INSTALLDIR)
	install -m 644 $(MODULES) $(INSTALLDIR)


clean:


uninstall:
	-rm -rf $(INSTALLDIR)
