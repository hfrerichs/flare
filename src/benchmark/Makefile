SUBDIRS = mesh_spacing
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))


.PHONY: all
all: $(SUBDIRS)


.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@


.PHONY: clean
clean: $(SUBCLEAN)


.PHONY: $(SUBCLEAN)
$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
