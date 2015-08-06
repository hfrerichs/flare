DIRS = plot_grid deposition_grid


.PHONY: all install clean $(DIRS)
all: $(DIRS)


install: $(DIRS)


clean: $(DIRS)


$(DIRS):
	cd $@; $(MAKE) $(MAKECMDGOALS)


addons:
	$(MAKE) all
addons_install:
	$(MAKE) install
addons_clean:
	$(MAKE) clean
