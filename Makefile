include include.mk


.PHONY: all
all:
	cd src; $(MAKE)


.PHONY: install
install:
	ln -s $(PWD)/bin/run_flare.sh $(BIN_DIR)
	ln -s $(PWD)/bin/$(PROGRAM) $(BIN_DIR)
	ln -s $(PWD)/bin/$(PROGRAM_DEBUG) $(BIN_DIR)


.PHONY: clean
clean:
	cd src; $(MAKE) clean

addon_targets = addons addons_install addons_clean
.PHONY: $(addon_targets)
$(addon_targets):
	cd src; $(MAKE) $(MAKECMDGOALS)
