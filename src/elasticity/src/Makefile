SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

CC = mpicc
CXX = mpicxx -fpermissive
LINK = mpicc

MKL_inc	= $(MKLROOT)/include

LOC_inc = ../../../local_include
UTS_inc = ../../utils/include
MAT_inc = ../../material/include
 	     	  
INCs = -I../include -I${LOC_inc} -I${UTS_inc} -I${MAT_inc} -I${MKL_inc}

CXXFLAGS = -Wall -std=c++14 -O3 -g
# for BGQ machine
#CFLAGS = -std=c99 -O3 -g -DARCH_BGQ

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
