include ../include.mk

SCRIPT = flare
TARGETS = flare-info flare-plot flare-equi2d flare-equi3d flare-geometry


.PHONY: all install clean uninstall $(TARGETS)
all: $(TARGETS)


$(TARGETS):
	sed "s|PYTHON|$(PYTHON)|" $@.py > $@


install:
	install -d $(BINDIR)
	install $(PROGRAM) $(BINDIR)
	install $(SCRIPT) $(BINDIR)
	install $(TARGETS) $(BINDIR)


clean:


uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(BINDIR)/$(SCRIPT)
	-rm -f $(foreach target,$(TARGETS),$(BINDIR)/$(target))
