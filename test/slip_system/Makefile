SRC = $(shell ls *.c)
OBJs = $(SRC:.c=.o)
PROG = test_slip_system

#########################################################
# set path
#########################################################
DIR_TOP = $(PWD)/../..
DIR_uts = $(DIR_TOP)/utils
DIR_mat = $(DIR_TOP)/material
DIR_els = $(DIR_TOP)/elasticity
DIR_cpm = $(DIR_TOP)/crystal_plasticity
DIR_loc = $(DIR_TOP)/local_include

CNSTVM_lib = -L$(DIR_TOP)/lib -lConstitutiveModel

CC = mpicc
CXX = mpicxx
LINK = mpicc

AR       = ar
ARFLAGS = -rcs

MKL_inc	= $(MKLROOT)/include
MKL_lib		= -Wl,--start-group \
 	     	  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_sequential.a \
 	     	  $(MKLROOT)/lib/intel64/libmkl_core.a \
 	     	  -Wl,--end-group -lpthread 	     	  

INCs = -I$(DIR_uts)/include \
	-I$(DIR_mat)/include \
	-I$(DIR_els)/include \
	-I$(DIR_cpm)/include \
	-I$(DIR_loc) \
	-I${MKL_inc}

LIBs = $(CNSTVM_lib) $(MKL_lib)
CFLAGS = -Wall -std=c99 -O3 -g

all: $(OBJs)
	$(LINK) $(CFLAGS) $(OBJs) $(LIBs) -o $(PROG)  
	
%.o:%.c
	$(CC) $(CFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) $(PROG) *.o
	
