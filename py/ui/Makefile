include ../../include.mk

MODULES = __init__.py grid.py


.PHONY: all install clean uninstall
all:


install:
	install -d $(PYTHONDIR)/flare/ui
	install -m 644 $(MODULES) $(PYTHONDIR)/flare/ui


clean:


uninstall:
	-rm -rf $(PYTHONDIR)/flare/ui
