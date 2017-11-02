include ../../include.mk

TARGETS = main_ui.py

MODULES = __init__.py flareui.py


.PHONY: all install clean uninstall
all: $(TARGETS)


%.py: %.ui
	pyside-uic $< > $@


install:
	install -d $(PYTHONDIR)/flare/gui
	install -m 644 $(TARGETS) $(MODULES) $(PYTHONDIR)/flare/gui


clean:
	-rm -f $(TARGETS)


uninstall:
	-rm -rf $(PYTHONDIR)/flare/ui
