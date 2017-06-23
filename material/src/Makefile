SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

CXX = mpicxx -fpermissive

MKL_inc	= $(MKLROOT)/include
 	     	  
INCs = -I../include -I../../local_include -I../../utils/include -m64 -I${MKL_inc}
CXXFLAGS = -Wall -std=c++14 -O3 -g
# for BGQ machine
#CFLAGS = -std=c99 -O3 -g -DARCH_BGQ

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
