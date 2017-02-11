SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)
PROG = test_J2_plasticity

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_uts = $(DIR_TOP)/src/utils
DIR_mat = $(DIR_TOP)/src/material
DIR_els = $(DIR_TOP)/src/elasticity
DIR_cpm = $(DIR_TOP)/src/crystal_plasticity
DIR_loc = $(DIR_TOP)/local_include
DIR_j2p = $(DIR_TOP)/src/J2_plasticity

CNSTVM_lib = -L$(DIR_TOP)/lib -lConstitutiveModel

CC = mpicc
CXX = mpicxx -fpermissive
LINK = mpicc

AR       = ar
ARFLAGS = -rcs

MKL_inc	= $(MKLROOT)/include
MKL_lib		= -Wl,--start-group \
 	     	  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_sequential.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_core.a \
 	     	  -Wl,--end-group -lpthread 	     	  

INCs = -I$(DIR_uts)/include \
	-I$(DIR_mat)/include \
	-I$(DIR_els)/include \
	-I$(DIR_cpm)/include \
	-I$(DIR_loc) \
	-I$(DIR_j2p)/include \
	-I${MKL_inc}

LIBs = $(CNSTVM_lib) $(MKL_lib)
CXXFLAGS = -Wall -std=c++14 -O3 -g

all: $(OBJs)
	$(LINK) $(CXXFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
