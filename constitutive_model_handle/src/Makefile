SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

CXX = mpicxx

#MKL_inc	= $(MKLROOT)/include

#LOC_inc = ../../include
#UTS_inc = ../utils/include
#MAT_inc = ../material/include
#ELS_inc = ../elasticity/include
#DAG_inc = ../damage/include
#PLS_inc = ../crystal_plasticity/include

#INCs = -I../include -I${LOC_inc} -I${UTS_inc} -I${MAT_inc} -I${ELS_inc} -I${DAG_inc} -I${PLS_inc} -I${MKL_inc}
INCs = -I../include

CXXFLAGS = -Wall -std=c++14 -O3
# for BGQ machine
#CFLAGS = -std=c99 -O3 -DARCH_BGQ

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
