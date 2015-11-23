# define compiler
CC = g++
# define compile time flags
CFLAGS = -g3 -Wall -fopenmp
# define any directories over than /usr/include/ containing header files
CINCLUDES = -I/usr/local/include -I . -I/usr/local/include/yaml-cpp
# define libraries paths over than /usr/include
LDFLAGS = -lfftw3 -lfftw3_omp -lm -lyaml-cpp
# define cpp source files (currently not using)
# SRCS = diffuseFT1D.cpp fftw++.cc (currently not using)
# make objects with macro (currently not using)
# OBJS = $(SRCS:.cpp = .o)
# OBJs += $(SRCS:.cc = .o)

OBJS = HardRodMain.o fftw++.o spGrid.o tGrid.o Propagator.o Ml3dPrint2Screen.o 
OBJS +=	ParseParams.o MayerFncHardRod.o

MAIN = diffuse

.PHONY: all

all: $(MAIN)

$(MAIN): $(OBJS) 
	$(CC) $(CINCLUDES) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) 

.PHONY: clean
clean:
	rm *.o

# $@: target $^: all dependenies $<: first dependency
