include include.mk

alldirs = src pyplot bin

.PHONY: all
all:
	cd src; $(MAKE)

.PHONY: debug
debug:
	cd src; $(MAKE) debug

.PHONY: install
install: $(alldirs)
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
uninstall: $(alldirs)
	-rm -rf $(DATADIR)/DIII-D/mockup_1
	-rm -rf $(DATADIR)/DIII-D/mockup_2
	-rm -rf $(DATADIR)/ITER/mockup_1
	-rm -rf $(DATADIR)/NSTX/mockup_1
	-rm $(DATADIR)/DIII-D/vessel_mockup.txt
	-rm $(DATADIR)/ITER/vessel_mockup.txt
	-rm $(DATADIR)/NSTX/vessel_mockup.txt
	-rm -d $(DATADIR)


.PHONY: $(alldirs)
$(alldirs):
	cd $@; $(MAKE) $(MAKECMDGOALS)


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd src; $(MAKE) $(MAKECMDGOALS)
