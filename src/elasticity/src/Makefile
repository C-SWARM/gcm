SRC = $(shell ls *.c)

OBJs = $(SRC:.c=.o)

CC = mpicc
CXX = mpicxx
LINK = mpicc

MKL_inc	= $(MKLROOT)/include

LOC_inc = ../../../local_include
UTS_inc = ../../utils/include
MAT_inc = ../../material/include
 	     	  
INCs = -I../include -I${LOC_inc} -I${UTS_inc} -I${MAT_inc} -I${MKL_inc}

CFLAGS = -Wall -std=c99 -O3 -g
# for BGQ machine
#CFLAGS = -std=c99 -O3 -g -DARCH_BGQ

all: $(OBJs) $(OLIB)
	
%.o:%.c
	$(CC) $(CFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
