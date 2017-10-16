include ../../include.mk


OBJECTS =  \
               magnetic_axis.o\
               reconstruct.o\
               bfield_component.o\
               equilibrium_format.o\
               m3dc1.o\
               geqdsk.o\
               divamhd.o\
               amhd.o\
               sonnet.o\
               equilibrium.o\
               mfc_polygon.o\
               polygones.o\
               splineB.o\
               interpolateB.o\
               libHINT.o\
               HINT.o\
               toroidal_harmonics.o\
               bfield.o

INC_DIR = -I .. -I ../$(EXTERNAL_DIR) -I ../$(CORE_DIR)


.PHONY: all
all: $(OBJECTS)

m3dc1.o: m3dc1.F90
	$(FC) -c $(INC_DIR) $(M3DC1_FLAG) $(M3DC1_INC) $<

reconstruct.o: reconstruct.F90
	$(FC) -c $(INC_DIR) $(EMC3_FLAG) -I ../$(EMC3_LINK_DIR) $<

%.o: %.f
	$(FC) -c $(INC_DIR) $<

%.o: %.f90
	$(FC) -c $(INC_DIR) $<

%.o: %.F90
	$(FC) -c $(INC_DIR) $<

.PHONY: clean
clean:
	rm -rf *.o
	rm -rf *.mod
