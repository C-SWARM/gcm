include ../../make.inc

SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

MKL_inc	= $(MKLROOT)/include
 	     	  
INCs = -I../include -I../../local_include -I../../utils/include -I${MKL_inc}

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
