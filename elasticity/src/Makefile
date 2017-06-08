include ../../make.inc

SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

LOC_inc = ../../local_include
DAT_inc = ../../data_structure
UTS_inc = ../../utils/include
MAT_inc = ../../material/include
 	     	  
INCs = -I../include -I${LOC_inc} -I${DAT_inc} -I${UTS_inc} -I${MAT_inc} -I${MKL_inc}

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
