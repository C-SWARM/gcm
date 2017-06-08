include ../../make.inc

SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

CXX = mpicxx

INCs = -I../include

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
