include ../../make.inc

SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)
PROG = test_data_structure

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_dat = $(DIR_TOP)/data_structure

INCs = -I$(DIR_dat)/include -I$(MKL_inc)

LIBs = $(MKL_lib)

all: $(OBJs)
	$(LINK) $(CXXFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
