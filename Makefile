OLIB = libConstitutiveModel.a

#########################################################
# compiler
#########################################################
CC = mpicc
CXX = mpicxx
LINK = mpicc

AR       = ar
ARFLAGS = -rcs

#########################################################
# get objs
#########################################################

SRC_uts = $(shell ls utils/src/*.cc)
OBJs_uts = $(SRC_uts:.cc=.o)

SRC_mat = $(shell ls material/src/*.cc)
OBJs_mat = $(SRC_mat:.cc=.o)

SRC_els = $(shell ls elasticity/src/*.cc)
OBJs_els = $(SRC_els:.cc=.o)

SRC_cpm = $(shell ls crystal_plasticity/src/*.cc)
OBJs_cpm = $(SRC_cpm:.cc=.o)

SRC_dam = $(shell ls damage/src/*.cc)
OBJs_dam = $(SRC_dam:.cc=.o)

SRC_cmh = $(shell ls constitutive_model_handle/src/*.cc)
OBJs_cmh = $(SRC_cmh:.cc=.o)

SRC_j2p = $(shell ls J2_plasticity/src/*.cc)
OBJs_j2p = $(SRC_j2p:.cc=.o)

SRC_pvp = $(shell ls poro_viscoplasticity/src/*.cc)
OBJs_pvp = $(SRC_pvp:.cc=.o)

OBJs = $(OBJs_uts) $(OBJs_mat) $(OBJs_els) $(OBJs_cpm) $(OBJs_dam) $(OBJs_cmh) $(OBJs_j2p) $(OBJs_pvp)

all: SRC


SRC:
	@echo $(OBJs)
	cd utils/src; make; cd ../..;
	cd material/src; make; cd ../..;
	cd elasticity/src; make; cd ../..;
	cd crystal_plasticity/src; make; cd ../..;
	cd damage/src; make; cd ../..;
	cd constitutive_model_handle/src; make; cd ../..;
	cd J2_plasticity/src; make; cd ../..;
	cd poro_viscoplasticity/src; make; cd ../..;

	$(AR) $(ARFLAGS) $(OLIB) $(OBJs)
	mv $(OLIB) lib	

clean:
	cd utils/src; make clean; cd ../..;
	cd material/src; make clean; cd ../..;
	cd elasticity/src; make clean; cd ../..;
	cd crystal_plasticity/src; make clean; cd ../..;
	cd damage/src; make clean; cd ../..;
	cd constitutive_model_handle/src; make clean; cd ../..;
	cd J2_plasticity/src; make clean; cd ../..;
	cd poro_viscoplasticity/src; make clean; cd ../..;
	cd lib; rm *.a; cd ..;
