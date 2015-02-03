include ../../include.mk


OBJECTS =  \
               magnetic_axis.o\
               reconstruct.o\
               geqdsk.o\
               divamhd.o\
               equilibrium.o\
               polygones.o\
               m3dc1.o\
               splineB.o\
               interpolateB.o\
               bfield.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(EMC3_LINK_DIR) -I ../$(CORE_DIR)


.PHONY: all
all: $(OBJECTS)

m3dc1.o: m3dc1.F90
	$(FC) -c $(INC_DIR) $(M3DC1_FLAG) $(M3DC1_INC) $<

%.o: %.f*
	$(FC) -c $(INC_DIR) $<

%.o: %.F*
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
