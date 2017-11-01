include ../include.mk

MODULES = __init__.py plot_lc.py


.PHONY: all install clean uninstall
all:


install:
	install -d $(PYTHONDIR)/flare/pyplot
	install -m 644 $(MODULES) $(PYTHONDIR)/flare/pyplot


clean:


uninstall:
	-rm -f $(PYTHONDIR)/flare/pyplot
