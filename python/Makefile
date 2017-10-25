include ../include.mk

alldirs = flare


.PHONY: all install clean uninstall $(alldirs)
all: $(alldirs)


install: $(alldirs)
	install -d $(BINDIR)
	install $(FLAREUI) $(BINDIR)


clean: $(alldirs)


uninstall:
	-rm -f $(BINDIR)/$(FLAREUI)


$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)
