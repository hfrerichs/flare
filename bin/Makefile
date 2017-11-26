include ../include.mk

SCRIPT = flare
TARGETS = flare-info flare-plot flare-equi2d flare-grid


.PHONY: all install clean uninstall $(TARGETS)
all:


$(BINDIR):
	install -d $(BINDIR)


$(TARGETS):
	install $@.py $(BINDIR)/$@


install: $(BINDIR) $(TARGETS)
	install $(PROGRAM) $(BINDIR)
	install $(SCRIPT) $(BINDIR)


clean:


uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(BINDIR)/$(SCRIPT)
	-rm -f $(foreach target,$(TARGETS),$(BINDIR)/$(target))
