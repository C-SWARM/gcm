SRC = $(shell ls *.c)

OBJs = $(SRC:.c=.o)

CC = mpicc
CXX = mpicxx
LINK = mpicc

MKL_inc	= $(MKLROOT)/include
 	     	  
INCs = -I../include -I../../local_include -I../../utils/include -I${MKL_inc}
CFLAGS = -Wall -std=c99 -O3 -g
# for BGQ machine
#CFLAGS = -std=c99 -O3 -g -DARCH_BGQ

all: $(OBJs) $(OLIB)
	
%.o:%.c
	$(CC) $(CFLAGS) $(INCs) -c $<

clean :
	rm -rf $(OLIB) *.o
	
