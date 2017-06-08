include ../../make.inc

SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)
PROG = test_crystal_plasticity_with_thermal

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_uts = $(DIR_TOP)/utils
DIR_mat = $(DIR_TOP)/material
DIR_els = $(DIR_TOP)/elasticity
DIR_cpm = $(DIR_TOP)/crystal_plasticity
DIR_loc = $(DIR_TOP)/local_include

CNSTVM_lib = -L$(DIR_TOP)/lib -lConstitutiveModel

INCs = -I$(DIR_uts)/include \
	-I$(DIR_mat)/include \
	-I$(DIR_els)/include \
	-I$(DIR_cpm)/include \
	-I$(DIR_loc) \
	-I${MKL_inc}

LIBs = $(CNSTVM_lib) $(MKL_lib)

all: $(OBJs)
	$(LINK) $(CXXFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o *.txt
	
