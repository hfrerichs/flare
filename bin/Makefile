include ../include.mk

SCRIPT = flare


.PHONY: all install clean uninstall
all:


install:
	install -d $(BINDIR)
	install $(PROGRAM) $(BINDIR)
	install $(SCRIPT) $(BINDIR)
	install $(RUN_FLARE) $(BINDIR)
	install $(FLAREUI) $(BINDIR)


clean:


uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(BINDIR)/$(SCRIPT)
	-rm -f $(BINDIR)/$(RUN_FLARE)
	-rm -f $(BINDIR)/$(FLAREUI)
