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
	test -d $(DATA_DIR) || mkdir $(DATA_DIR)


.PHONY: clean
clean:
	cd src; $(MAKE) clean


.PHONY: uninstall
uninstall:
	-rm -f $(BIN_DIR)/run_flare.sh
	-rm -f $(BIN_DIR)/$(PROGRAM)
	-rm -f $(BIN_DIR)/$(PROGRAM_DEBUG)
	-rm -d $(BIN_DIR)
	-rm -d $(DATA_DIR)


addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd src; $(MAKE) $(MAKECMDGOALS)
