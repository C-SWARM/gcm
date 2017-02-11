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

SRC_uts = $(shell ls src/utils/*.cc)
OBJs_uts = $(SRC_uts:.cc=.o)

SRC_mat = $(shell ls src/material/*.cc)
OBJs_mat = $(SRC_mat:.cc=.o)

SRC_els = $(shell ls src/elasticity/*.cc)
OBJs_els = $(SRC_els:.cc=.o)

SRC_cpm = $(shell ls src/crystal_plasticity/*.cc)
OBJs_cpm = $(SRC_cpm:.cc=.o)

SRC_dam = $(shell ls src/damage/*.cc)
OBJs_dam = $(SRC_dam:.cc=.o)

SRC_cmh = $(shell ls src/constitutive_model_handle/*.cc)
OBJs_cmh = $(SRC_cmh:.cc=.o)

SRC_j2p = $(shell ls src/J2_plasticity/*.cc)
OBJs_j2p = $(SRC_j2p:.cc=.o)

OBJs = $(OBJs_uts) $(OBJs_mat) $(OBJs_els) $(OBJs_cpm) $(OBJs_dam) $(OBJs_cmh) $(OBJs_j2p)

all: SRC


SRC:
	@echo $(OBJs)
	cd src/utils; make; cd ../..;
	cd src/material; make; cd ../..;
	cd src/elasticity; make; cd ../..;
	cd src/crystal_plasticity; make; cd ../..;
	cd src/damage; make; cd ../..;
	cd src/constitutive_model_handle; make; cd ../..;
	cd src/J2_plasticity; make; cd ../..;

	$(AR) $(ARFLAGS) $(OLIB) $(OBJs)
	mv $(OLIB) lib	

clean:
	cd src/utils; make clean; cd ../..;
	cd src/material; make clean; cd ../..;
	cd src/elasticity; make clean; cd ../..;
	cd src/crystal_plasticity; make clean; cd ../..;
	cd src/damage; make clean; cd ../..;
	cd src/constitutive_model_handle; make clean; cd ../..;
	cd src/J2_plasticity; make clean; cd ../..;
	cd lib; rm *.a; cd ..;
