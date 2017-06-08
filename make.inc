
# compile
CXX = mpicxx -fpermissive
LINK = mpicxx
AR       = ar
ARFLAGS = -rcs

# compile options
CXXFLAGS = -Wall -std=c++14 -O3
# for BGQ machine
#CFLAGS = -std=c99 -O3 -g -DARCH_BGQ

# libraries
MKL_inc = $(MKLROOT)/include

MKL_lib = $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
 $(MKLROOT)/lib/intel64/libmkl_sequential.a \
 $(MKLROOT)/lib/intel64/libmkl_core.a -lpthread -ldl

