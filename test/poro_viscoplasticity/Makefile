SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)
PROG = poro_viscoplasticity

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_uts = $(DIR_TOP)/utils
DIR_mat = $(DIR_TOP)/material
DIR_els = $(DIR_TOP)/elasticity
DIR_cpm = $(DIR_TOP)/crystal_plasticity
DIR_loc = $(DIR_TOP)/local_include
DIR_pvp = $(DIR_TOP)/poro_viscoplasticity

CNSTVM_lib = -L$(DIR_TOP)/lib -lConstitutiveModel

CXX = mpicxx
LINK = mpicxx

MKL_inc	= $(MKLROOT)/include
MKL_lib		= -Wl,--start-group \
 	     	  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_sequential.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_core.a \
 	     	  -Wl,--end-group -lpthread -lm -ldl -lstdc++

INCs = -I$(DIR_pvp)/include \
	-I$(DIR_mat)/include \
	-I${MKL_inc}

LIBs = $(CNSTVM_lib) $(MKL_lib)
CXXFLAGS = -Wall -std=c++14 -O3 -g

all: $(OBJs)
	$(LINK) $(CXXFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o out
	
