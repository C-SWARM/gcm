include ../../make.inc

SRC = $(shell ls *.cc)

OBJs = $(SRC:.cc=.o)

MKL_inc	= $(MKLROOT)/include

LOC_inc = ../../local_include
UTS_inc = ../../utils/include
MAT_inc = ../../material/include
ELS_inc = ../../elasticity/include
 	     	  
INCs = -I../include -I${LOC_inc} -I${UTS_inc} -I${MAT_inc} -I${ELS_inc}  -I${MKL_inc}

all: $(OBJs) $(OLIB)
	
%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
