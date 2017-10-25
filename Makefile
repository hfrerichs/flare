include include.mk


.PHONY: all
all:
	cd src; $(MAKE)

.PHONY: debug
debug:
	cd src; $(MAKE) debug

.PHONY: install
install:
	test -d $(BINDIR) || install -d $(BINDIR)
	test -d $(LIBDIR) || install -d $(LIBDIR)

	install bin/run_flare.sh $(BINDIR)
	install bin/$(PROGRAM) $(BINDIR)
	install src/$(FLARELIB) $(LIBDIR)
	mkdir -p $(DATADIR)/DIII-D
	mkdir -p $(DATADIR)/ITER
	mkdir -p $(DATADIR)/NSTX
	cp templates/Database/DIII-D/vessel_mockup.txt $(DATADIR)/DIII-D
	cp -r templates/Database/DIII-D/mockup_1 $(DATADIR)/DIII-D
	cp -r templates/Database/DIII-D/mockup_2 $(DATADIR)/DIII-D
	cp templates/Database/ITER/vessel_mockup.txt $(DATADIR)/ITER
	cp -r templates/Database/ITER/mockup_1 $(DATADIR)/ITER
	cp templates/Database/NSTX/vessel_mockup.txt $(DATADIR)/NSTX
	cp -r templates/Database/NSTX/mockup_1 $(DATADIR)/NSTX


.PHONY: clean
clean:
	cd src; $(MAKE) clean


.PHONY: uninstall
uninstall:
	-rm -f $(BINDIR)/run_flare.sh
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -d $(BINDIR)
	-rm -f $(LIBDIR)/$(FLARELIB)
	-rm -d $(LIBDIR)
	-rm -rf $(DATADIR)/DIII-D/mockup_1
	-rm -rf $(DATADIR)/DIII-D/mockup_2
	-rm -rf $(DATADIR)/ITER/mockup_1
	-rm -rf $(DATADIR)/NSTX/mockup_1
	-rm $(DATADIR)/DIII-D/vessel_mockup.txt
	-rm $(DATADIR)/ITER/vessel_mockup.txt
	-rm $(DATADIR)/NSTX/vessel_mockup.txt
	-rm -d $(DATADIR)


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd src; $(MAKE) $(MAKECMDGOALS)
