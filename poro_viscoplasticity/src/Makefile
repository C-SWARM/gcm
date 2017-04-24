SRC = $(shell ls *.cc)
OBJs = $(SRC:.cc=.o)

#########################################################
# set compiler and path
#########################################################

CXX = g++ 

CXXFLAGS = -Wall -std=c++14 -O3 -g

MAT_inc = ../../material/include

INCs = -I../include -I$(MKLROOT)/include -I$(MAT_inc)

all: $(OBJs)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf *.o

