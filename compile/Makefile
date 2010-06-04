
# Default Makefile

BASE = ..
NAME = cm3d

# common definitions...
RM = rm -f
CC = cc
F77 = gfortran
CPP = cc -E
CP = cp
AR = ar
MV = mv

# Compiler flages, paths for include files and libraries

#DEFINES = -DLINUX -DBLAS -DFORTRANUNDERSCORE
#DEFINES = -DPRESSURE 
DEFINES = -DLINUX 

CFLAGS = -g -Wall $(DEFINES) 
FFLAGS =   

#OPT = -O
OPT =
LIBS = -lm
#LIBS = -lm -lg2c -L/usr/lib/libblas -L/usr/lib/liblapack

LD   = $(CC)  
ARFLAGS = ruv
RANLIB = ranlib

include Makefile.cm3d

