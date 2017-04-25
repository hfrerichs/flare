include include.mk


.PHONY: all
all:
	cd src; $(MAKE)

.PHONY: debug
debug:
	cd src; $(MAKE) debug

.PHONY: install
install:
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	ln -s $(PWD)/bin/run_flare.sh $(BIN_DIR)
	ln -s $(PWD)/bin/$(PROGRAM) $(BIN_DIR)
	ln -s $(PWD)/bin/$(PROGRAM_DEBUG) $(BIN_DIR)
	mkdir -p $(DATA_DIR)/DIII-D
	mkdir -p $(DATA_DIR)/ITER
	mkdir -p $(DATA_DIR)/NSTX
	cp templates/Database/DIII-D/vessel_mockup.txt $(DATA_DIR)/DIII-D
	cp -r templates/Database/DIII-D/mockup_1 $(DATA_DIR)/DIII-D
	cp -r templates/Database/DIII-D/mockup_2 $(DATA_DIR)/DIII-D
	cp templates/Database/ITER/vessel_mockup.txt $(DATA_DIR)/ITER
	cp -r templates/Database/ITER/mockup_1 $(DATA_DIR)/ITER
	cp templates/Database/NSTX/vessel_mockup.txt $(DATA_DIR)/NSTX
	cp -r templates/Database/NSTX/mockup_1 $(DATA_DIR)/NSTX


.PHONY: clean
clean:
	cd src; $(MAKE) clean


.PHONY: uninstall
uninstall:
	-rm -f $(BIN_DIR)/run_flare.sh
	-rm -f $(BIN_DIR)/$(PROGRAM)
	-rm -f $(BIN_DIR)/$(PROGRAM_DEBUG)
	-rm -d $(BIN_DIR)
	-rm -rf $(DATA_DIR)/DIII-D/mockup_1
	-rm -rf $(DATA_DIR)/DIII-D/mockup_2
	-rm -rf $(DATA_DIR)/ITER/mockup_1
	-rm -rf $(DATA_DIR)/NSTX/mockup_1
	-rm $(DATA_DIR)/DIII-D/vessel_mockup.txt
	-rm $(DATA_DIR)/ITER/vessel_mockup.txt
	-rm $(DATA_DIR)/NSTX/vessel_mockup.txt
	-rm -d $(DATA_DIR)


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd src; $(MAKE) $(MAKECMDGOALS)
