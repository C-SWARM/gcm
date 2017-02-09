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

SRC_uts = $(shell ls src/utils/src/*.c)
OBJs_uts = $(SRC_uts:.c=.o)

SRC_mat = $(shell ls src/material/src/*.c)
OBJs_mat = $(SRC_mat:.c=.o)

SRC_els = $(shell ls src/elasticity/src/*.c)
OBJs_els = $(SRC_els:.c=.o)

SRC_cpm = $(shell ls src/crystal_plasticity/src/*.c)
OBJs_cpm = $(SRC_cpm:.c=.o)

SRC_dam = $(shell ls src/damage/src/*.c)
OBJs_dam = $(SRC_dam:.c=.o)

SRC_cmh = $(shell ls src/constitutive_model_handle/src/*.c)
OBJs_cmh = $(SRC_cmh:.c=.o)

SRC_j2p = $(shell ls src/J2_plasticity/src/*.c)
OBJs_j2p = $(SRC_j2p:.c=.o)

OBJs = $(OBJs_uts) $(OBJs_mat) $(OBJs_els) $(OBJs_cpm) $(OBJs_dam) $(OBJs_cmh) $(OBJs_j2p)

all: SRC


SRC:
	@echo $(OBJs)
	cd src/utils/src; make; cd ../../..;
	cd src/material/src; make; cd ../../..;
	cd src/elasticity/src; make; cd ../../..;
	cd src/crystal_plasticity/src; make; cd ../../..;
	cd src/damage/src; make; cd ../../..;
	cd src/constitutive_model_handle/src; make; cd ../../..;
	cd src/J2_plasticity/src; make; cd ../../..;

	$(AR) $(ARFLAGS) $(OLIB) $(OBJs)
	mv $(OLIB) lib	

clean:
	cd src/utils/src; make clean; cd ../../..;
	cd src/material/src; make clean; cd ../../..;
	cd src/elasticity/src; make clean; cd ../../..;
	cd src/crystal_plasticity/src; make clean; cd ../../..;
	cd src/damage/src; make clean; cd ../../..;
	cd src/constitutive_model_handle/src; make clean; cd ../../..;
	cd src/J2_plasticity/src; make clean; cd ../../..;
	cd lib; rm *.a; cd ..;
