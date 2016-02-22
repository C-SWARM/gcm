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

SRC_uts = $(shell ls utils/src/*.c)
OBJs_uts = $(SRC_uts:.c=.o)

SRC_mat = $(shell ls material/src/*.c)
OBJs_mat = $(SRC_mat:.c=.o)

SRC_els = $(shell ls elasticity/src/*.c)
OBJs_els = $(SRC_els:.c=.o)

SRC_cpm = $(shell ls crystal_plasticity/src/*.c)
OBJs_cpm = $(SRC_cpm:.c=.o)

SRC_dam = $(shell ls damage/src/*.c)
OBJs_dam = $(SRC_dam:.c=.o)

SRC_cmh = $(shell ls constitutive_model_handle/src/*.c)
OBJs_cmh = $(SRC_cmh:.c=.o)

OBJs = $(OBJs_uts) $(OBJs_mat) $(OBJs_els) $(OBJs_cpm) $(OBJs_dam) $(OBJs_cmh)

all: SRC


SRC:
	@echo $(OBJs)
	cd utils/src; make; cd ../..;
	cd material/src; make; cd ../..;
	cd elasticity/src; make; cd ../..;
	cd crystal_plasticity/src; make; cd ../..;
	cd damage/src; make; cd ../..;
	cd constitutive_model_handle/src; make; cd ../..;	
	$(AR) $(ARFLAGS) $(OLIB) $(OBJs)
	mv $(OLIB) lib	

clean:
	cd utils/src; make clean; cd ../..;
	cd material/src; make clean; cd ../..;
	cd elasticity/src; make clean; cd ../..;
	cd crystal_plasticity/src; make clean; cd ../..;
	cd damage/src; make clean; cd ../..;
	cd constitutive_model_handle/src; make clean; cd ../..;
	cd lib; rm *.a; cd ..;
