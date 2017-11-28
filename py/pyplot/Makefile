include ../../include.mk

MODULES = __init__.py plot_lc.py plot_1d.py plot_2d.py


.PHONY: all install clean uninstall
all:


install:
	install -d $(PYTHONDIR)/flare/pyplot
	install -m 644 $(MODULES) $(PYTHONDIR)/flare/pyplot


clean:


uninstall:
	-rm -rf $(PYTHONDIR)/flare/pyplot
