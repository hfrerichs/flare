include ../../../include.mk


INC_DIR = -I ../..

.PHONY: all
all: $(EMC3_OBJ)


%.o: %.f
	$(FC) -c $(INC_DIR) $<

PHYS_CONST.f:
	ln -s $(EMC3_SRC_DIR)/../USER/PHYS_CONST.f

SURFACE_PL.f:
	ln -s $(EMC3_SRC_DIR)/SURFACE_PL.f

MAPPING.f:
	ln -s $(EMC3_SRC_DIR)/MAPPING.f

check.f:
	sed -n '/CHECK_GRID_POINTS/,/CHECK_GRID_POINTS/p'   $(EMC3_SRC_DIR)/check.f > check.f
	sed -n '/PLOT_CROSS_SECTION/,/PLOT_CROSS_SECTION/p' $(EMC3_SRC_DIR)/check.f >> check.f

ibm_iface.f:
	ln -s $(EMC3_SRC_DIR)/../SERV/ibm_iface.f

random.f:
	ln -s $(EMC3_SRC_DIR)/../SERV/random.f

service.f:
	ln -s $(EMC3_SRC_DIR)/../SERV/service.f

sf_def_user.f:
	ln -s $(EMC3_SRC_DIR)/../USER/sf_def_user.f

sf_jump.f:
	ln -s $(EMC3_SRC_DIR)/sf_jump.f

.PHONY: clean
clean:
	rm -rf *.o *.mod
	for file in *f; do if [ -L $$file ]; then rm $$file; fi; done
	rm check.f
