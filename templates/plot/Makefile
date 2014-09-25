OBJECTS =          plot_lc.pro\
                   plot_data.pro


.PHONY: all
all: $(OBJECTS)


plot_lc.pro: grid.src plot_lc.src
	cat grid.src > plot_lc.pro
	cat plot_lc.src >> plot_lc.pro


plot_data.pro: grid.src plot_data.src
	cat grid.src > plot_data.pro
	cat plot_data.src >> plot_data.pro


.PHONY: clean
clean:
	rm -rf plot_lc.pro plot_data.pro
