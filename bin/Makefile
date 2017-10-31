include ../include.mk


.PHONY: all install clean uninstall
all:


install:
	install -d $(BINDIR)
	install $(PROGRAM) $(BINDIR)
	install $(RUN_FLARE) $(BINDIR)
	install $(FLAREUI) $(BINDIR)


clean:


uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(BINDIR)/$(RUN_FLARE)
	-rm -f $(BINDIR)/$(FLAREUI)
