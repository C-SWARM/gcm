SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)
PROG = test_data_structure

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_dat = $(DIR_TOP)/data_structure

CXX = mpicxx
LINK = mpicxx

MKL_inc	= $(MKLROOT)/include
MKL_lib	= $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -lpthread -ldl 

INCs = -I$(DIR_dat)/include -I$(MKL_inc)

LIBs = $(MKL_lib)
CXXFLAGS = -Wall -std=c++14 -O3

all: $(OBJs)
	$(LINK) $(CXXFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
