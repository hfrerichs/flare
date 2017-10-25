include ../../../include.mk

TARGETS = main_ui.py

MODULES = __init__.py flareui.py

.PHONY: all install clean uninstall
all:


install:
	install -d $(PYTHONDIR)/flare/ui
	install -m 644 $(TARGETS) $(MODULES) $(PYTHONDIR)/flare/ui


clean:
	-rm -f $(TARGETS)


uninstall:
	-rm -rf $(PYTHONDIR)/flare/ui
