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

CELL_GEO.f:
	ln -s $(EMC3_SRC_DIR)/CELL_GEO.f

cell_def.f:
	ln -s $(EMC3_SRC_DIR)/cell_def.f

cell_user.f:
	ln -s $(EMC3_SRC_DIR)/../USER/cell_user.f

check.f:
	ln -s $(EMC3_SRC_DIR)/check.f

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
