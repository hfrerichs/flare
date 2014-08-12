include ../../include.mk


OBJECTS =  \
               trace_bline.o\
               sample_bfield.o\
               poincare_plot.o\
               connection_length.o\
               get_equi_info_2D.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(BFIELD_DIR) -I ../$(CORE_DIR) -I ../$(GEOMETRY_DIR)


.PHONY: all
all: $(OBJECTS)

%.o: %.f*
	$(FC) -c $(INC_DIR) $<

%.o: %.F*
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
