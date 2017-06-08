include ../../make.inc

SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

INCs = -I./include -I../include -I${MKL_inc}

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
